% load('mycamp')
% load position.mat
load('BCICIV_1\electrodesBCICIV1.mat')
load HeadModel.mat
rel = rand(1,59);
sel = logical(ones(1,59));
rel = rel./max(abs(rel));
MyTopo_fun3(rel,sel,M1.xy,M1.lab,[min(rel) max(rel)],0,0)
caxis([-1,1])
% colormap(gca,cmap);
axis square
axis off
