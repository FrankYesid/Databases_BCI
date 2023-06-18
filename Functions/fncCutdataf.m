function Xdr = fncCutdataf(X,tini,tfin,fs,filterbands) 

%% recortar (MI segment)
ini = tini*fs;%0.1*fs;
fin = tfin*fs;
% Xdr = cell(1,length(Xd_));
%for i = 1:length(Xd_)
    %Xf_= fcnfiltband(Xd_{i},fs,[8,30],5);
    %X = Xd_;
    for k = 1:length(X)
        [samples,~] = size(X{k});
        lags = 1:samples;  
        %pos = find(lags<fin+1 & lags>ini);
        Xtmp{1,k} = X{k}(lags<fin+1 & lags>ini,:);  
    end
    
    Xdr = fcnfiltband(Xtmp,fs,filterbands,5);
%end 