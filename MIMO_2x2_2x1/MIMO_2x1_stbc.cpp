/*
---------------------------------------------MIMO 2X1 ALAMOUTI STBC, CON MODULACIÓN QAM------------------------------------------------
La simulación consiste en lo siguiente:
    -Generador de Bits aleatorios
    -Modulador QAM
    -Codificación y decodificación de Alamouti con ruido AWGN con diferente SNR y canal Rayleigh
    -Demodulador QAM
    -Medición de BER

Los resultados que entrada la simulación son:
    -Un vector con los valores  SNR en dB utilizados
    -Un vector con los resultados de BER obtenidos
    -Un archivo MIMO_2X1_result_file.it con los valores obtenidos para realizar la graficación en
      MATLAB u Octave
*/


#include <itpp/itcomm.h>
using namespace itpp;
//Estas lineas son para usar el cout y endl de la libreria Itpp
using std::cout;
using std::endl;


int main()
{
//int
int M = 256;         //Tamaño de la constelación (M = 16,64,256)
int k = log2(M);    //Bits por símbolo
int N = 1000000;    //Número de bits a simular
int i, j, z;
int y;
int qamLength;


//Double
double Ec,Eb;

//Vectores
vec SNRdB, BER;
vec EbN0dB, EbN0, N0;

//Vectores binarios
bvec in_bits, dec_bits;

//Vectores complejos
cvec qamSymbols, r0n, r1n;
cvec s0, s1, h0;
cvec h1, r0, r1;
cvec h0neg, reciv;

//Matrices complejas
cmat H, R, Hh;
cmat recR0R1;


//Classes
QAM qam(M);                         //Se define el modulador M-QAM
BERC berc;                          //Clase de contador de BER
AWGN_Channel awgn_channel(0.1);     //Clase de canal AWGN
it_file ff;                         //Clase para fuardar resultados en archivo


//Parámetros del ruido AWGN
Ec = 1.0;                       //La energía transmitida por símbolo QAM es 1
Eb = Ec / k;                    //La energía transmitida por bit es Ec/tamaño de palabra
SNRdB = linspace(0,20,21);      //Vector de SNR en dB
EbN0dB = SNRdB - 10*log10(k);   //Se convierte el parámetro SNR a N0, para poder utilizar
EbN0 = inv_dB(EbN0dB);          //la clase AWGN predefinida en la libreria
N0 = Eb * pow(EbN0, -1.0);


//Aleatorización de los generadores de valores rand
RNG_randomize();

//Se generan N-bits aleatorios:
in_bits = randb(N);

//Se modula el vector bits y se guardan los símbolos en el vector qamSymbols
qam.modulate_bits(in_bits, qamSymbols);
qamLength = qamSymbols.length();        //Tamaño del vector de símbolos

//Dividimos el vector de símbolos en s0 y s1, un ejemplo a continuación:
//qamSymbols = [a b c d e f];
//s0 = [a c e]
//s1 = [b d f]
s0 = zeros_c(1,qamLength/2);
s1 = zeros_c(1,qamLength/2);
z = 0;
for (j=0; j < qamLength; j = j+2){
    if (j == 0){
        s0(j) = qamSymbols(j);
        s1(j) = qamSymbols(j+1);
    }else {
        s0(z) = qamSymbols(j);
        s1(z) = qamSymbols(j+1);
    }
    z++;
}


//Caracterización del canal h0 y h1, canales Rayleigh
std::complex<double> mycomplex (0,1);
h0 = 1/sqrt(2)*(randn(1,qamLength/2)) + mycomplex * (randn(1,qamLength/2));
h1 = 1/sqrt(2)*(randn(1,qamLength/2)) + mycomplex * (randn(1,qamLength/2));


//Señal recibida sin ruido
r0 = elem_mult(h0,s0) + elem_mult(h1,s1);               //r0 = h0*s0 + h1*s1
h0neg = -1*h0;
r1 = elem_mult(h0neg,conj(s1)) + elem_mult(h1,conj(s0));//r1 = -h0*conj(s1) + h1*conj(s0)


BER = zeros(1,SNRdB.length());                      //Vector para guardar los resultados de BER totales

//Ciclo de codificación de Alamouti con ruido, decodificación de Alamouti,
//demodulación QAM y cálculo de BER
for (i=0; i < SNRdB.length(); i++) {

    //Señal recibida con ruido
    awgn_channel.set_noise(N0(i));
    r0n = awgn_channel(r0);             //r0 = h0*s0 + h1*s1 + n0
    r1n = awgn_channel(r1);             //r1 = -h0*conj(s1) + h1*conj(s0) + n1

    H = zeros_c(2,2);                   //Vector vacío para guardar la matriz de canal
    R = zeros_c(2,1);                   //Vector vacío para guardaar la matriz de señal recibida
    reciv = zeros_c(1,(r0.length()+r1.length()));   //Vector para almacenar los símbolos recividos

    //------------------Decodificación de Alamouti----------------------------------------------
    y = 0;
    for (k = 0; k < h0.length(); k++) {
        //Formamos la matriz de canal H = [h0 h1;h1* -h0*]
        H(0,0) = h0(k);
        H(0,1) = h1(k);
        H(1,0) = conj(h1(k));
        H(1,1) = -1*conj(h0(k));
        //Calculamos la matriz pseudoinversa
        Hh = conj(H);
        Hh = Hh.transpose();
        H = inv(Hh*H)*Hh;       //Matriz pseudo inversa, equivalente a H = pinv(H) en MATLAB

        R(0,0) = r0n(k);
        R(1,0) = conj(r1n(k));  //Vector de señal recivida R = [r0; conj(r1)]

        recR0R1 = H*R;          //Decodificación de símbolos
        recR0R1 = recR0R1.transpose();
        //Símbolos decodificados guardados en la variable reciv
        reciv(y) = recR0R1(0,0);
        reciv(y+1) = recR0R1(0,1);
        y = y + 2;
    }
    //----------------Fin de decodificación de Alamouti-------------------------------------------

    //Se demodula el vector qamSymbols y se guardan esos bits en dec_bits
    qam.demodulate_bits(reciv, dec_bits);
    //cout<<"Dec bits = "<< dec_bits <<endl;

    //Cálculo de BER
    berc.clear();
    berc.count(in_bits, dec_bits);

    //Se guardan los datos en el vector de resultados BER
    BER(i) = berc.get_errorrate();

}
//Se muestran los resultados totales
cout<<"------------------------Total--------------------------------"<<endl;
cout<<"SNR dB ="<< SNRdB <<endl;
cout<<"BER ="<< BER <<endl;


//Se almacenan los resultados de la simulación en un archivo .it:
ff.open("MIMO_2X1_result_file.it");
ff << Name("SNR_2x1") << SNRdB;
ff << Name("BER_2x1") << BER;
ff.close();
}
