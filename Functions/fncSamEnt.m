function [SampEn] = fncSamEnt(series,ul1,up1,wind)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function for different channels
%
%

ws = ul1:wind:up1;
SEtemp = zeros(1,length(ws));
parfor a = 1:length(ws)
    temp = series(ws(a):ws(a)+wind-1);
    SEtemp(a) = SampEn2(2,0.2*std(temp),temp);
end
SampEn =  SEtemp;



% fprintf('SampEnt Channel - %d. Time: %f\n',ch,time)
