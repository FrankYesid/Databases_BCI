function [] = fncFFT(X1,X2,Fs)
%fourier
L = length(X1);

Y = fft(X1);
P2 = abs(Y/L);
P1_ = P2(1:L/2+1);
P1_(2:end-1) = 2*P1_(2:end-1);

Y = fft(X2);
P2 = abs(Y/L);
P2_ = P2(1:L/2+1);
P2_(2:end-1) = 2*P2_(2:end-1);


f = Fs*(0:(L/2))/L;
loglog(f,P1_/max(max([P1_;P2_])));
hold on
loglog(f,P1_/max(max([P1_;P2_])));
    
end