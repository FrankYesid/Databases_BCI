function [Xout] = fncrotate(W,varargin)
% Rota los datos y los Hankeliza
%
X = varargin;
Xout = cell(1,nargin-1);
for i = 1: (nargin-1)
    Xtmp = X{i};
    ntr = size(Xtmp,2);
    [~,covtmp] = fncCov(Xtmp,0);
    Xrot = cell(1,ntr);
    for tr = 1:ntr
        Xrot{tr} = Xtmp{tr}./sqrt(trace(covtmp(:,:,tr)));
        Xrot{tr} = W*Xtmp{tr}';
    end
    Xout{i} = Xrot;
end

% 
% 
% [~,covtr] = fncCov(Xtrain);
% 
% ntr = size(Xtrain,2);
% nts = size(Xtest,2);
% Xrtrain = cell(1,ntr);
% Xrtest = cell(1,nts);
% 
% for tr = 1:ntr
% %     [trace(covtr(:,:,tr))]
%     Xtrain{tr} = Xtrain{tr}./sqrt(trace(covtr(:,:,tr)));
%     Xrtrain{tr} = W*Xtrain{tr}';
% end
% 
% clear covtr
% [~,covtr] = fncCov(Xtest);
% for tr = 1:nts
%     Xtest{tr} = Xtest{tr}./sqrt(trace(covtr(:,:,tr)));
%     Xrtest{tr} = W*Xtest{tr}';
% end

