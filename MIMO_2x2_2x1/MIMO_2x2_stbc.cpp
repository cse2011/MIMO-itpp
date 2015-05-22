/*
---------------------------------------------MIMO 2X2 ALAMOUTI STBC, CON MODULACIÓN QAM------------------------------------------------
La simulación consiste en lo siguiente:
    -Generador de Bits aleatorios
    -Modulador QAM
    -Codificación y decodificación de Alamouti con ruido AWGN con diferente SNR y canal Rayleigh
    -Demodulador QAM
    -Medición de BER

Los resultados que entrada la simulación son:
    -Un vector con los valores  SNR en dB utilizados
    -Un vector con los resultados de BER obtenidos
    -Un archivo MIMO_2X2_result_file.it con los valores obtenidos para realizar la graficación en
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
int N = 3000000;    //Número de bits a simular
int i, j, z;
int y, a, hLength;
int qamLength;


//Double
double Ec,Eb;

//Vectores
vec SNRdB, BER;
vec EbN0dB, EbN0, N0;

//Vectores binarios
bvec in_bits, dec_bits;

//Vectores complejos
cvec qamSymbols, r00n, r11n;
cvec r01n, r10n;
cvec s0, s1, h00;
cvec h11, r0, r1;
cvec h10, h01;
cvec h00neg, reciv;
cvec r00;
cvec r01;
cvec r10;
cvec r11;

//Matrices complejas
cmat H, R, Hh;
cmat recR0R1, re, H2;
cmat S;


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
N0 = Eb * pow(EbN0, -1.0);  ;


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


//Caracterización del canal h00 y h11, canales Rayleigh
std::complex<double> mycomplex (0,1);
h00 = 1/sqrt(2)*(randn(1,qamLength/2)) + mycomplex * (randn(1,qamLength/2));
h10 = 1/sqrt(2)*(randn(1,qamLength/2)) + mycomplex * (randn(1,qamLength/2));
h11 = 1/sqrt(2)*(randn(1,qamLength/2)) + mycomplex * (randn(1,qamLength/2));
h01 = 1/sqrt(2)*(randn(1,qamLength/2)) + mycomplex * (randn(1,qamLength/2));

hLength = h00.length();

//Señal recibida sin ruido

r00 = zeros_c(1,hLength);
r01 = zeros_c(1,hLength);
r10 = zeros_c(1,hLength);
r11 = zeros_c(1,hLength);

H2 = zeros_c(4,2);
S = zeros_c(2,1);

for (a=0; a < hLength; a++){

    H2(0,0) = h00(a);
    H2(0,1) = h01(a);
    H2(1,0) = h10(a);
    H2(1,1) = h11(a);
    H2(2,0) = conj(h01(a));
    H2(2,1) = -1*conj(h00(a));
    H2(3,0) = conj(h11(a));
    H2(3,1) = -1*conj(h10(a));

    S(0,0) = s0(a);
    S(1,0) = s1(a);

    re = H2*S;

    r00(a) = re(0,0);
    r01(a) = re(1,0);
    r10(a) = re(2,0);
    r11(a) = re(3,0);
}



BER = zeros(1,SNRdB.length());                      //Vector para guardar los resultados de BER totales

//Ciclo de codificación de Alamouti con ruido, decodificación de Alamouti,
//demodulación QAM y cálculo de BER
for (i=0; i < SNRdB.length(); i++) {

    //Señal recibida con ruido
    awgn_channel.set_noise(N0(i));
    r00n = awgn_channel(r00);             //r0 = h00*s0 + h11*s1 + n0
    r01n = awgn_channel(r01);             //r0 = h00*s0 + h11*s1 + n0
    r10n = awgn_channel(r10);             //r0 = h00*s0 + h11*s1 + n0
    r11n = awgn_channel(r11);             //r1 = -h00*conj(s1) + h11*conj(s0) + n1

    H = zeros_c(4,2);                   //Vector vacío para guardar la matriz de canal
    R = zeros_c(4,1);                   //Vector vacío para guardaar la matriz de señal recibida
    reciv = zeros_c(1,(r00.length()+r11.length()));   //Vector para almacenar los símbolos recividos

    //------------------Decodificación de Alamouti----------------------------------------------
    y = 0;
    for (k = 0; k < h00.length(); k++) {
        H = zeros_c(4,2);
        //Formamos la matriz de canal H = [h00 h11;h11* -h00*]
        H(0,0) = h00(k);
        H(0,1) = h01(k);
        H(1,0) = h10(k);
        H(1,1) = h11(k);
        H(2,0) = conj(h01(k));
        H(2,1) = -1*conj(h00(k));
        H(3,0) = conj(h11(k));
        H(3,1) = -1*conj(h10(k));

        //Calculamos la matriz pseudoinversa
        Hh = conj(H);
        Hh = Hh.transpose();
        H = inv(Hh*H)*Hh;       //Matriz pseudo inversa, equivalente a H = pinv(H) en MATLAB

        R(0,0) = r00n(k);
        R(1,0) = r01n(k);  //Vector de señal recivida R = [r0; conj(r1)]
        R(2,0) = r10n(k);
        R(3,0) = r11n(k);


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
ff.open("MIMO_2X2_result_file.it");
ff << Name("SNR_2x2") << SNRdB;
ff << Name("BER_2x2") << BER;
ff.close();
}
