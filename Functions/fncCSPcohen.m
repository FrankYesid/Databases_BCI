function [W,eigenvalues] = fncCSPcohen(Xd,Xdr,mode,flag)
%Input:
%Xd : X data
%Xr : X reference

%Output:
%W: the learnt CSP filters (a [Nc*Nc] matrix with the filters as rows)

[covXd,~]= fncCov(Xd,mode,0);               % Data covariance
[covXr,~] = fncCov(Xdr,mode,1);             % Reference covariance


% covXd = covXd./max(max(covXd));
% covXr = covXr./max(max(covXr));

% covXd = (log(abs(covXd)));
% covXr = (log(abs(covXr)));

% figure('name','covXd'),imagesc(covXd),title('CovXd'),colorbar
% figure('name','covXr'),imagesc(covXr),title('CovXr'),colorbar

%whitening transform of total covariance matrix
[Ut,Dt] = eig(covXr); %caution: the eigenvalues are initially in increasing order
% [Ut,Dt] = eig(covXd+covXr);
% eig(covXd+covXr);
eigenvalues = diag(Dt);
[eigenvalues,egIndex] = sort(eigenvalues, 'descend');
Ut = Ut(:,egIndex);
P = diag(sqrt(1./(eigenvalues))) * Ut';

%transforming covariance matrix of first class using P
transformedCov1 =  P * (covXd) * P';

%EVD of the transformed covariance matrix
[U1,D1] = eig(transformedCov1);
eigenvalues = diag(D1);
[eigenvalues,egIndex] = sort(eigenvalues, 'descend');
% U = U1;

if flag ==0
%     U1 = U1(:, egIndex);
    U1 = U1(:, egIndex([end]));
elseif flag ==1
    U1 = U1(:, egIndex([1 end]));
elseif flag ==3
    U1 = U1(:, egIndex([1]));
end
W = U1' * P;

%W = W;