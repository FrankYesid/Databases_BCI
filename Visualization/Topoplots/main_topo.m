% clear all; close all; clc
function main_topo(rel)
load('BCICIV_2a/electrodesBCICIV2a.mat')
load('HeadModel.mat')   % model of the head.
Fil = 3; Col = 3;       % topoplot en subplots para cada uno de los sujetos.
t_e = 40;               % tamaño visual de los electrodos seleccionados.
M1.xy = elec_pos(:,1:2);% posicion de los canales.
M1.lab = Channels;      % nombre de los canales.
% MM = M1.lab(1:11,1);
% MM2 = M1.lab(12:end,1);
% M1.lab = [MM2;MM];
% mask = reshape(rho,[17 22]) <= thresh; %

 %%
sel = 1:22;
%                 sel = find(mask(band,:)==1); % Orden de los canales segun una relevancia.
%                 tmp = W{fold,band}(caract,:); % rel = linspace(1,10,numel(M1.lab));     % para colocar mas importancia al electrodo.
% load('C:\Users\frany\Downloads\rel.mat');
% aa = rel_a.rel(1:11,1);
% bb = rel_a.rel(12:end,1);
% rel = [bb;aa];
% rel = ones(1,22);
% rel = rel/max(abs(rel));
pos = M1.xy;

label = M1.lab;
warning off
for i=1:2
    pos(:,i) = 0.9.*((pos(:,i)-min(pos(:,i)))/(range(pos(:,i)))-0.5);
end
xc = HeadModel(1,:);
yc = HeadModel(2,:);
% figure
hold on
% Topoplot
x = pos(:,1);
y = pos(:,2);
tmp = [x,y,x*0]*rotz(2);
x = tmp(:,1);
y = tmp(:,2);
pos(:,1) = x;
pos(:,2) = y;
GS = 500;
xi = linspace(min(x)-0.05, max(x)+0.05, GS);       % x-axis for interpolation (row vector)
yi = linspace(min(y)-0.05, max(y)+0.05, GS);       % y-axis for interpolation (row vector)
[Xi, Yi, Zi] = griddata([x' xc], [y' yc], [rel(:);zeros(numel(xc),1)], xi', yi,'v4'); % interpolate the topographic data
%% Creating data mask
[TH R] = cart2pol(Xi,Yi);
Zi(R>0.5) = NaN;
deltax = xi(2)-xi(1); % length of grid entry
deltay = yi(2)-yi(1); % length of grid entry

h = surf(Xi-deltax/2, Yi-deltay/2+0.004, zeros(size(Zi)), Zi,'EdgeColor', 'none', 'FaceColor', 'flat');hold on
caxis([0 1])
colormap('parula')
% shading interp
%                 plot(0.99.*xc,0.99.*yc,'Color',[130,130,130]/255,'LineWidth',6)
tam = 1;
plot(tam*xc,tam*yc,'color',[156,156,156]./255,'LineWidth',2)
scatter(0.9.*x(sel),0.9.*y(sel),t_e,'b','filled')
% view([-90 90])
% text(pos(sel,1)*0.9-0.02,pos(sel,2)*0.9-0.04,label(sel),'Interpreter',...
%     'latex','ButtonDownFcn',{@lineCallback,database},'Color','white',...
%     'FontSize',10);
axis off
hold off
axis square



