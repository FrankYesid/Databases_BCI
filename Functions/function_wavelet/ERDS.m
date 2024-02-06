function [ERD] = ERDS(x,t1,t2)
% x dimensiones (trial,filter,ch,tiempo)
t = 0:1/250:7-(1/250);
temp1 = abs(t - t1);
min1 = min(temp1);
temp2 = abs(t - t2);
min2 = min(temp2);
ul = find(temp1 == min1);
up = find(temp2 == min2);

% 
r_nc = squeeze(mean(x(:,:,:,ul:up),4));

r_c = squeeze(mean(r_nc,1));

m_c = squeeze(mean(x,1));

% ERD 
ERD = bsxfun(@times,m_c,1./r_c) - 1;
