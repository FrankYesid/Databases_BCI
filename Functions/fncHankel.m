function [Xrtrain,Xrtest] = fncHankel(Xtrain,Xtest,L)
% Rota los datos y los Hankeliza
%

ntr = size(Xtrain,2);
nts = size(Xtest,2);
Xrtrain = cell(1,ntr);
Xrtest = cell(1,nts);

for tr = 1:ntr
    tmp = hankel(Xtrain{tr});
    Xrtrain{tr} = tmp(:,1:L);
%     Xrtrain{tr} = Xrtrain{tr}(:,1:250);
end

for tr = 1:nts
    tmp = hankel(Xtest{tr});
    Xrtest{tr} = tmp(:,1:L);
%     Xrtest{tr} = Xrtest{tr}(:,1:250);
end
