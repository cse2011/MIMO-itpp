figure(1); clf;
itload('MIMO_2X1_alamouti_qam_result_file.it')
itload('MIMO_2X1_alamouti_qpsk_result_file.it')
h1 = semilogy(SNR_QAM,BER_QAM,'*-r'); hold on;
h2 = semilogy(SNR_QPSK,BER_QPSK,'+-b');
set(gca,'fontname','times','fontsize',16);
xlabel('SNR [dB]');
ylabel('BER')
title('Comparación del rendimiento de la BER con modulación QAM y QPSK');
legend([h1 h2],'16 QAM','QPSK');
grid on;

%p = path
%path(p,'/usr/share/itpp')