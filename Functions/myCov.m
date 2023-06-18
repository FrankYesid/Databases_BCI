function C=myCov(X)

m = mean(X,1);
temp = bsxfun(@minus,X,m);
C = temp'*temp/(size(temp,1)-1);


