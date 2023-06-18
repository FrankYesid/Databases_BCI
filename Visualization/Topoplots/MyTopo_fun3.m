function MyTopo_fun3(Y,sel,pos,label,lims,cur,ticks)
% Y vector de pesos 1xN
% pos vector de posiciones Nx2
% label celda de etiquetas Nx1
% lims vector limites del colorbar [min max]
% cur booleana habilita curvas de nivel 1,0 
% ticks habilida plotear etiquetas boolena 1,0

%Y = Y/max(Y);


for i=1:2
    pos(:,i) = 0.9*((pos(:,i)-min(pos(:,i)))/range(pos(:,i))-0.5);
end
%figure
%HeadModel
load HeadModel
xc = HeadModel(1,:);
yc = HeadModel(2,:);

plot(xc,yc,'k','LineWidth',3)
hold on
%% Topoplot
GS = 100;
x = pos(:,1);
y = pos(:,2);
xi         = linspace(min(x)-0.05, max(x)+0.05, GS);       % x-axis for interpolation (row vector)
yi         = linspace(min(y)-0.05, max(y)+0.05, GS);       % y-axis for interpolation (row vector)
% [Xi, Yi, Zi] = griddata(x', y, Y, xi', yi,'v4'); % interpolate the topographic data
[Xi, Yi, Zi] = griddata([x' xc], [y' yc], [Y(:);zeros(numel(xc),1)], xi', yi,'v4'); % interpolate the topographic data
%% Creating data mask
[TH R] = cart2pol(Xi,Yi);
Zi(R>0.5) = NaN;
deltax = xi(2)-xi(1); % length of grid entry
deltay = yi(2)-yi(1); % length of grid entry
h = surf(Xi-deltax/2, Yi-deltay/2, zeros(size(Zi)), Zi,'EdgeColor', 'none', 'FaceColor', 'flat');hold on

% Color = repmat([0 0 1],numel(Y(sel)>0),1);
% Color(Y(sel)>0,:) = repmat([1 1 0],sum(Y(sel)>0),1);

scatter(x(sel),y(sel),sel(sel)*400,[0 0 0],'filled')
hold on
if cur == 1
    contour(Xi-deltax/2, Yi-deltay/2,Zi,'k')
    hold on
end

plot(xc,yc,'k','LineWidth',5)

if ticks == 1
    ii = ismember(label,'')
    plot(x(not(ii)),y(not(ii)),'*')
    set(gca,'XTick',[],'YTick',[]);
    text(x,y+0.02,label);
end

% aa = colorbar('ytick',0:5);
%set(aa,'YTick',0:5)
caxis(lims)

function [xunit yunit] = circle(x,y,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
