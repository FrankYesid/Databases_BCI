function [SampEn] = Main_SE(series,ul1,up1,t3)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function for different channels
%
%
[N_trial,~,channel] = size(series);
ws = up1-ul1+1;
SampEn = zeros(N_trial,ws,channel);

for ch = 1 : channel
    tic
    for trial = 1:N_trial
%          tic
        dat = squeeze(series(trial,:,ch));
        parfor a = 1:ws
            ti = a + ul1-1;
            temp = dat(ti:ti+t3);
            SEtemp(a) = SampEn2(2,0.45*std(temp),temp);
        end
        SampEn(trial,:,ch) =  SEtemp;
        
%          fprintf('Trial - %d. Time: %f\n',trial,time)
    end
    time = toc;
    fprintf('SampEnt Channel - %d. Time: %f\n',ch,time)
end