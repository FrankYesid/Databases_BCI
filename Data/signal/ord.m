function signal = ord(X) 
% X: trials.

for tr = 1:numel(X)
    N_channel = size(X{tr},2);
    for ch = 1:N_channel
      signal{ch,1}(:,tr) = X{tr}(:,ch);  
    end
end
