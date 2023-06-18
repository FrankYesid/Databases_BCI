%% Main para graficar la entropias

% entropy
entropy = rand(2,10,22);
%% parametros
% limites de la amplitud de la señal 
limi = [0,1];
% labels del eje X
labelx = num2str([0,1]);
% labels del eje Y
labely = num2str([0,1]);
sel = 1:22;
limits_y = [0 1];
limits_x = [1 size(entropy,2)];
tam_label = 6;
wi = 0.1;
hi = 0.07;
% tamaño de la grafica
set(gcf,'position',[667   528   404   420])
x1 = 0.02;     y1 = 0.06;     w = 0.95;      h = 0.92;
subplot('position',[x1,y1,w,h]);
labels = 0; % activo los numero de los ejes
main_topo_plots(entropy,limits_x,limits_y,sel,tam_label,wi,hi,labels)

% save
saveas(gcf,['head_s'],'epsc')