function [Xfreq] = fcnfiltband2(X,fs,freq,n)
%% FUNCTION rhythms
%-------------------------------------------------------------------------
% Inputs:
% Xf = A_filterButter(Xcell,fs)
% X{trial}(samples x chan): Data to analysis
% 
% fs: Sample frequency
% n : filter order
% freq : vector con las bandas
% 
%-------------------------------------------------------------------------
% Outputs:
% Xf: Filtered data 
% Xf{trial}(samples x chan)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2017 Signal Processing and Recognition Group
% L.F. Velasquez-Martinez 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%  filter order
if nargin < 4  
   n = 5; % by default 
end

%Filter desing
wn_freq = [freq(1)/(0.5*fs) freq(2)/(0.5*fs)]; %
[b_freq,a_freq] = butter(n,wn_freq);  % filtrado

% show filter response
%freqz(b_alpha,a_alpha,512,fs)
%title('n = 5 Butterworth band pass Filter')
tr = length(X);
Xfreq = cell(tr,1);

for k = 1: tr % channel
    [s,c] = size(X{k});
    if c>s X{k} = X{k}'; [s,c] = size(X{k}); end
    xpru = cell2mat(X{k});
    for ch = 1 : s
        x = xpru(ch,:);
        Xfreq{k,1}{ch,1} = filtfilt(b_freq,a_freq,x); %filter
    end
end
