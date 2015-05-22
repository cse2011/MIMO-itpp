clear
clc

itload('MIMO_2X1_result_file.it')
itload('MIMO_2X2_result_file.it')
h1 = semilogy(SNR_2x1,BER_2x1,'*-r'); hold on
h2 = semilogy(SNR_2x2,BER_2x2,'+-b');
title("MIMO 2x1 STBC 16-QAM");
xlabel("SNR dB"); 
ylabel("BER");
legend([h1 h2],'MIMO 2X1 16-QAM','MIMO 2X2 16-QAM');
grid on;

