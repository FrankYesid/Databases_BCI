function [tt] = fncTest(A,B)
[N_channel,~] = size(A);
parfor chan = 1:N_channel
    [~,p] = ttest2(A(chan,:),B(chan,:));
    %     tl(1,1) = h;
    %     tl(1,2) = p;
    tt(chan,:) = p;
end


% function [tt] = tes2(A,B)
% [trials,~,channels] = size(A);
% tt = zeros(trials,channels);
% for tr = 1:trials
%     for ch = 1:channels
%         [~,p] = ttest2(A(tr,:,ch),B(tr,:,ch));
% %         tl(1,1) = h;
% %         tl(1,2) = p;
%         tt(tr,ch) = p;
%     end
% end