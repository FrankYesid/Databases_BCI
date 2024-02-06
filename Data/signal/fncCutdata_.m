
function Xdr = fncCutdata_(Xd_,tini,tfin,fs) 

%% recortar (MI segment)
ini = tini*fs;%0.1*fs;
fin = tfin*fs;
Xdr = cell(1,length(Xd_));
% for i = 1:length(Xd_)
    Xf_ = Xd_;
    for k = 1:length(Xf_)
        [samples,~] = size(Xf_{k});
        lags = 1:samples;  
        pos = find(lags<fin+1 & lags>ini);
        Xdr{k} = Xf_{k}(pos,:);  
    end
% end