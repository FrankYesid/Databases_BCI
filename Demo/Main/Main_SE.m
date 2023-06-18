%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function for different channels
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [SampEn] = Main_SE(series,ul1,up1,t3)
[ch,f,T] = size(series);
series = reshape(series,[ch*f,T]);
SampEn = nan(ch*f,up1-ul1+1);
parfor i = 1:ch*f
    a=1;
    temSE = zeros(1,up1-ul1+1);
    for ti = ul1:up1
        temSE(1,a) = SampEn2(2,0.45*std(series(i,ti:ti+t3)),(series(i,ti:ti+t3)));
        a = a +1;
    end
     SampEn(i,:) = temSE;
end
SampEn = reshape(SampEn,[ch f up1-ul1+1]);
% [ch,f,~] = size(series);  
% tm = length(ul1:up1); 
% SampEn = zeros(ch,f,tm);
% tic
% for chan = 1:ch
%     
%     for freq = 1:f 
%         dat = squeeze(series(chan,freq,:));
%         %a=1;
%         for ti = ul1:up1
%             %            x = SampEn2(2,0.45*std(dat(ti:ti+t3)),dat(ti:ti+t3));
%             %            flag =  isinf(x);
%             %            if flag == 1
%             %                break
%             %            end
%             temSE(ti-100) = SampEn2(2,0.45*std(dat(ti:ti+t3)),dat(ti:ti+t3));
%             %a = a +1;
%         end
%         SampEn(chan,freq,:) = temSE;
%     end
% %     fprintf(['Channel: ' num2str(chan) num2str(toc) '\n'])
% end   
% fprintf(['Dur SamEntropy: ' num2str(toc) '\n'])