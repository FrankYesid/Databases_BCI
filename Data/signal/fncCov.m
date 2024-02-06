function [covXr,covtr] = fncCov(Xdr,mode,regularization)
%input 
% Data reference
% mode = 1 - samples, 0- channels
%Output:

%%
% [nchan,nsamples] = size(Xdr{1});
% if nchan>nsamples; [nchan,~] = size(Xdr{1}'); end
ntr = length(Xdr);
%covtr = zeros(nchan,nchan,ntr);
nan_ = 0;
for j = 1:ntr
    E = Xdr{j};
    if mode ~= 0; E = E'; end 
    covtr(:,:,j) = cov(E);
    covtr(:,:,j)= covtr(:,:,j)/trace(covtr(:,:,j));
    %regularization
    if regularization==1
        covtr(:,:,j) = covtr(:,:,j)+(0.001.*(diag(covtr(:,:,j))*diag(covtr(:,:,j))'));
        fl = sum(sum(isnan(covtr(:,:,j))));
        if fl==1
            nan_= nan_ +1;
        end
    end
end
covXr = mean(covtr(:,:,:),3);
