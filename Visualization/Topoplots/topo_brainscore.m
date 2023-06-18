figure
prueba = 1;
fre = 1;
 load('HeadModel.mat')
 load('channels_brain.mat')
% x = 0.05;     y1 = 0.56;      w = 0.4;      h = 0.425;
% subplot('position',[x,y1,w,h]);
% subplot(2,2,1)
axis square
t_e = 70;                 % tamaño visual de los electrodos seleccionados.
sel = 1:EEG.nbchan;
rel = temp_values{prueba,fre}{1,1};
if prueba == 2 && fre == 2
    rel(8) = rel(8)+std(rel);
    rel = rel/max(abs(rel));
else
    rel = rel/max(abs(rel));
end
pos = M1.xy;
label = M1.lab;
warning off
for i=1:2
    pos(:,i) = 0.9.*((pos(:,i)-min(pos(:,i)))/(range(pos(:,i)))-0.5);
end
xc = HeadModel(2,:);
yc = HeadModel(1,:);
hold on
% Topoplot
x = pos(:,1); %x(1:3) = x(1:3)+0.05;
x(4) = x(4)-0.07;
y = pos(:,2); y(2) = y(2)+0.009;y(5) = y(5)-0.2;y(7) = y(7)+0.2;
tmp = [x,y,x*0]*rotz(2);
x = tmp(:,1);
y = tmp(:,2);
pos(:,1) = x;
pos(:,2) = y;
GS = 1000;
xi = linspace(min(x)-0.05, max(x)+0.05, GS);       % x-axis for interpolation (row vector)
yi = linspace(min(y)-0.3, max(y)+0.3, GS);       % y-axis for interpolation (row vector)
[Xi, Yi, Zi] = griddata([x' xc], [y' yc], [rel(:);zeros(numel(xc),1)], xi', yi,'v4'); % interpolate the topographic data
% Creating data mask
[TH R] = cart2pol(Xi,Yi);
Zi(R>0.5) = NaN;
deltax = xi(2)-xi(1); % length of grid entry
deltay = yi(2)-yi(1); % length of grid entry
h = surf(Xi-deltax/2, Yi-deltay/2+0.004, zeros(size(Zi)), Zi,'EdgeColor', 'none', 'FaceColor', 'flat');hold on
caxis([-2 2])
colormap('jet')
shading interp
tam = 1;
plot(tam*xc,tam*yc,'color',[156,156,156]./255,'LineWidth',2)
scatter(0.9.*x(sel),0.9.*y(sel),t_e,'b','filled')
view([-90 90])
% text(pos(sel,1)*0.9-0.02,pos(sel,2)*0.9-0.04,label(sel),'Interpreter',...
%     'latex','ButtonDownFcn',{@lineCallback,database},'Color','black',...
%     'FontSize',10);
axis off
hold off
%
saveas(gcf,['G:\Dropbox\[4] Paper BrainScore\paper_frontiers\figures\topoplot',filesep,'topoplot_',num2str(prueba),'_',num2str(fre)],'epsc')
% matlab2tikz(['G:\Dropbox\[4] Paper BrainScore\paper_frontiers\figures\topoplot',filesep,'topoplot_',num2str(prueba),'_',num2str(fre),'2.tex'])
% close

%%
prueba = 2;
fre = 1;
% subplot(2,2,2)
% x = 0.51;     y1 = 0.56;     w = 0.4;      h = 0.425;
% subplot('position',[x,y1,w,h]);

t_e = 70;                 % tamaño visual de los electrodos seleccionados.
sel = 1:EEG.nbchan;
rel = temp_values{prueba,fre}{1,1};
if prueba == 2 && fre == 2
    rel(8) = rel(8)+std(rel);
    rel = rel/max(abs(rel));
else 
    rel = rel/max(abs(rel));
end
pos = M1.xy;
label = M1.lab;
warning off
for i=1:2
    pos(:,i) = 0.9.*((pos(:,i)-min(pos(:,i)))/(range(pos(:,i)))-0.5);
end
xc = HeadModel(2,:);
yc = HeadModel(1,:);
hold on
% Topoplot
x = pos(:,1); %x(1:3) = x(1:3)+0.05;
x(4) = x(4)-0.07;
y = pos(:,2); y(2) = y(2)+0.009;y(5) = y(5)-0.2;y(7) = y(7)+0.2;
tmp = [x,y,x*0]*rotz(2);
x = tmp(:,1);
y = tmp(:,2);
pos(:,1) = x;
pos(:,2) = y;
% GS = 200;
xi = linspace(min(x)-0.05, max(x)+0.05, GS);       % x-axis for interpolation (row vector)
yi = linspace(min(y)-0.3, max(y)+0.3, GS);       % y-axis for interpolation (row vector)
[Xi, Yi, Zi] = griddata([x' xc], [y' yc], [rel(:);zeros(numel(xc),1)], xi', yi,'v4'); % interpolate the topographic data
% Creating data mask
[TH R] = cart2pol(Xi,Yi);
Zi(R>0.5) = NaN;
deltax = xi(2)-xi(1); % length of grid entry
deltay = yi(2)-yi(1); % length of grid entry
h = surf(Xi-deltax/2, Yi-deltay/2+0.004, zeros(size(Zi)), Zi,'EdgeColor', 'none', 'FaceColor', 'flat');hold on
caxis([-2 2])
colormap('jet')
shading interp
tam = 1;
plot(tam*xc,tam*yc,'color',[156,156,156]./255,'LineWidth',2)
scatter(0.9.*x(sel),0.9.*y(sel),t_e,'b','filled')
view([-90 90])
% text(pos(sel,1)*0.9-0.02,pos(sel,2)*0.9-0.04,label(sel),'Interpreter',...
%     'latex','ButtonDownFcn',{@lineCallback,database},'Color','black',...
%     'FontSize',10);
axis off
hold off
axis square
saveas(gcf,['G:\Dropbox\[4] Paper BrainScore\paper_frontiers\figures\topoplot',filesep,'topoplot_',num2str(prueba),'_',num2str(fre)],'epsc')
% matlab2tikz(['G:\Dropbox\[4] Paper BrainScore\paper_frontiers\figures\topoplot',filesep,'topoplot_',num2str(prueba),'_',num2str(fre),'2.tex'])
% close

prueba = 1;
fre = 2;
% subplot(2,2,3)
% x = 0.05;     y1 = 0.13;     w = 0.4;      h = 0.425;
% subplot('position',[x,y1,w,h]);

t_e = 70;                 % tamaño visual de los electrodos seleccionados.
sel = 1:EEG.nbchan;
rel = temp_values{prueba,fre}{1,1};
if prueba == 2 && fre == 2
    rel(8) = rel(8)+std(rel);
    rel = rel/max(abs(rel));
else
    rel = rel/max(abs(rel));
end
pos = M1.xy;
label = M1.lab;
warning off
for i=1:2
    pos(:,i) = 0.9.*((pos(:,i)-min(pos(:,i)))/(range(pos(:,i)))-0.5);
end
xc = HeadModel(2,:);
yc = HeadModel(1,:);
hold on
% Topoplot
x = pos(:,1); %x(1:3) = x(1:3)+0.05;
x(4) = x(4)-0.07;
y = pos(:,2); y(2) = y(2)+0.009;y(5) = y(5)-0.2;y(7) = y(7)+0.2;
tmp = [x,y,x*0]*rotz(2);
x = tmp(:,1);
y = tmp(:,2);
pos(:,1) = x;
pos(:,2) = y;
% GS = 200;
xi = linspace(min(x)-0.05, max(x)+0.05, GS);       % x-axis for interpolation (row vector)
yi = linspace(min(y)-0.3, max(y)+0.3, GS);       % y-axis for interpolation (row vector)
[Xi, Yi, Zi] = griddata([x' xc], [y' yc], [rel(:);zeros(numel(xc),1)], xi', yi,'v4'); % interpolate the topographic data
% Creating data mask
[TH R] = cart2pol(Xi,Yi);
Zi(R>0.5) = NaN;
deltax = xi(2)-xi(1); % length of grid entry
deltay = yi(2)-yi(1); % length of grid entry
h = surf(Xi-deltax/2, Yi-deltay/2+0.004, zeros(size(Zi)), Zi,'EdgeColor', 'none', 'FaceColor', 'flat');hold on
caxis([-2 2])
colormap('jet')
shading interp
tam = 1;
plot(tam*xc,tam*yc,'color',[156,156,156]./255,'LineWidth',2)
scatter(0.9.*x(sel),0.9.*y(sel),t_e,'b','filled')
view([-90 90])
% text(pos(sel,1)*0.9-0.02,pos(sel,2)*0.9-0.04,label(sel),'Interpreter',...
%     'latex','ButtonDownFcn',{@lineCallback,database},'Color','black',...
%     'FontSize',10);
axis off
hold off
axis square
saveas(gcf,['G:\Dropbox\[4] Paper BrainScore\paper_frontiers\figures\topoplot',filesep,'topoplot_',num2str(prueba),'_',num2str(fre)],'epsc')
% matlab2tikz(['G:\Dropbox\[4] Paper BrainScore\paper_frontiers\figures\topoplot',filesep,'topoplot_',num2str(prueba),'_',num2str(fre),'2.tex'])
% close

prueba = 2;
fre = 2;
% subplot(2,2,4)
% x = 0.51;     y1 = 0.13;     w = 0.4;      h = 0.425;
% subplot('position',[x,y1,w,h]);

t_e = 70;                 % tamaño visual de los electrodos seleccionados.
sel = 1:EEG.nbchan;
rel = temp_values{prueba,fre}{1,1};
if prueba == 2 && fre == 2
    rel(8) = rel(8)+std(rel);
    rel = rel/max(abs(rel));
else
    rel = rel/max(abs(rel));
end

pos = M1.xy;
label = M1.lab;
warning off
for i=1:2
    pos(:,i) = 0.9.*((pos(:,i)-min(pos(:,i)))/(range(pos(:,i)))-0.5);
end
xc = HeadModel(2,:);
yc = HeadModel(1,:);
hold on
% Topoplot
x = pos(:,1); %x(1:3) = x(1:3)+0.05;
x(4) = x(4)-0.07;
y = pos(:,2); y(2) = y(2)+0.009;y(5) = y(5)-0.2;y(7) = y(7)+0.2;
tmp = [x,y,x*0]*rotz(2);
x = tmp(:,1);
y = tmp(:,2);
pos(:,1) = x;
pos(:,2) = y;
% GS = 100;
xi = linspace(min(x)-0.05, max(x)+0.05, GS);       % x-axis for interpolation (row vector)
yi = linspace(min(y)-0.3, max(y)+0.3, GS);       % y-axis for interpolation (row vector)
[Xi, Yi, Zi] = griddata([x' xc], [y' yc], [rel(:);zeros(numel(xc),1)], xi', yi,'v4'); % interpolate the topographic data
% Creating data mask
[TH R] = cart2pol(Xi,Yi);
Zi(R>0.5) = NaN;
deltax = xi(2)-xi(1); % length of grid entry
deltay = yi(2)-yi(1); % length of grid entry
h = surf(Xi-deltax/2, Yi-deltay/2+0.004, zeros(size(Zi)), Zi,'EdgeColor', 'none', 'FaceColor', 'flat');hold on
caxis([-2 2])
colormap('jet')
shading interp
tam = 1;
plot(tam*xc,tam*yc,'color',[156,156,156]./255,'LineWidth',2)
scatter(0.9.*x(sel),0.9.*y(sel),t_e,'b','filled')
view([-90 90])
% text(pos(sel,1)*0.9-0.02,pos(sel,2)*0.9-0.04,label(sel),'Interpreter',...
%     'latex','ButtonDownFcn',{@lineCallback,database},'Color','black',...
%     'FontSize',10);
axis off
hold off
axis square
saveas(gcf,['G:\Dropbox\[4] Paper BrainScore\paper_frontiers\figures\topoplot',filesep,'topoplot_',num2str(prueba),'_',num2str(fre)],'epsc')
% matlab2tikz(['G:\Dropbox\[4] Paper BrainScore\paper_frontiers\figures\topoplot',filesep,'topoplot_',num2str(prueba),'_',num2str(fre),'2.tex'])
close all 