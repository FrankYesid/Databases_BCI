function [SampEn] = Main_SE2(series,ul1,up1,t3)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function for different channels
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w = up1-ul1+1;
N_fil = numel(series);
for fr = 1:N_fil
    [~,N_channel] = size(series{fr});
    for chan = 1:N_channel
        dat = squeeze(series{fr}(:,chan));
%         tic 
        parfor a = 1:w
            ti = a+ul1-1;
            SampEn(fr,chan,a) = SampEn2(2,0.45*std(dat(ti:ti+t3-1)),dat(ti:ti+t3));
        end
%        time = toc;
%        fprintf('SampEn: Channel %d. Time: %f\n',chan,time)
    end
end