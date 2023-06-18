%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function for different channels
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Reny] = Main_RE(series,ul1,up1,t3)
[ch,f,T] = size(series);
series = reshape(series,[ch*f,T]);
Reny = nan(ch*f,up1-ul1+1);
parfor i = 1:ch*f
    a=1;
    tmp = zeros(1,up1-ul1+1);
    for ti = ul1:up1
        tmp(1,a) = alphaH(series(i,ti:ti+t3)',128,2);
        a = a +1;
    end
    Reny(i,:) = tmp;
    %fprintf(['Channel: ' num2str(chan) num2str(toc) '\n'])
end
Reny = reshape(Reny,[ch f up1-ul1+1]);