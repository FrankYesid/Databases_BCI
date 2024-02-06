function [P1,f] = fncFFT(X,Fs,normalized)
%fourier
L = length(X);
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;

if normalized ==1
    P1 = P1/max(P1);
end
    
% plot(f,P1/max(P1));
% title('Single-Sided Amplitude Spectrum')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')
end