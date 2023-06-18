% topoplot giga
% close all
% load('posiciones.mat')
% load('HeadModel.mat')   % model of the head.
% t_e = 70;               % tamaño visual de los electrodos seleccionados.
% sel = (1:64);           % Orden de los canales segun una relevancia.
% M1.xy = floor(M1.xy(:,[1 2]).*100); % posicion de los canales.
% M1.xy(28,1) = -105;     % posicion de los canales.
% % M1.xy(16,2) = M1.xy(53,2);
% % M1.xy(28,1) = -10;
% M1.lab = M1.lab;        % nombre de los canales.
% database = 3;
% % relevancia.
% rel = ones(1,numel(M1.lab));
% rel = rel/max(abs(rel));

clear;
clc;

load('BCICIV_2a\electrodesBCICIV2a.mat')
%load('BCICIV_2a\labels.mat')
%load('BCICIV_2a\layout.mat')
load('HeadModel.mat')   % model of the head.
Fil = 3; Col = 3;             % topoplot en subplots para cada uno de los sujetos.
t_e = 70;               % tamaño visual de los electrodos seleccionados.
M1.xy = elec_pos(:,1:2);% posicion de los canales.
M1.lab = Channels;      % nombre de los canales.
% mask = reshape(rho,[17 22]) <= thresh; %
sel = 1:22;%find(mask(band,:)==1); % Orden de los canales segun una relevancia.
% relevancia para los canales la variable es w.
% tmp = W{fold,band}(caract,:); % rel = linspace(1,10,numel(M1.lab));     % para colocar mas importancia al electrodo.
rel = ones(1,22);
% rel(sel) = ;
time_ = [0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5];
for sss = [2,8]
    load(['I:\results rhy\ERDs_',num2str(sss),'.mat'])
    for cl = 1:2
        for ch = 1:22
            er{cl}(ch,:,:) = r{1}{cl}.ERDS{ch}.erds;
        end 
    end
    for tt = 1:numel(time_);t_(tt) = find(r{1}{1}.t_plot == time_(tt)); end
    for cl = 1:2
        for tim = 1:length(time_)
            for fre = 1:17
                %                 axes('Position',[0.05 0.05 0.89 0.89]);
                figure;
                fig = gca;
                %                 fig.Parent.PaperPosition = [2.91 4.15 2.69 2.69];
                fig.Parent.OuterPosition = [1998 541 406 468];
                rel = squeeze(er{cl}(:,t_(tim),fre));
                pos = M1.xy;
                label = M1.lab;
                % tam - tamaño de la cabeza y posición de los electrodos.
                tam = 1; tam2 = 0.96; tam3 = 0.98;
                for i=1:2
                    pos(:,i) = tam2.*((pos(:,i)-min(pos(:,i)))/(range(pos(:,i)))-0.5);
                end
                xc = HeadModel(1,:);
                yc = HeadModel(2,:);
                load('G:\Dropbox\ERD\results_ERDfc_subjects\Mapas de colores\erdscolormap.mat')
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
                shading interp % interpola colores.                
%                 scatter(tam3.*x(sel),(tam3+0.04).*y(sel),t_e,'b','filled')
%                 text(pos(sel,1)*tam3-0.02,pos(sel,2)*tam3+0.02,label(sel),'Interpreter',...
%                     'latex','ButtonDownFcn',{@lineCallback,database},'Color','black',...
%                     'FontSize',10);
                plot(tam*xc,tam*yc+0.001,'k','LineWidth',3)
                colormap(erdcolormap)
                caxis([-1 1.5])
                axis off
                hold off
                axis square
                axis image;   
%                 fig = gca;
%                 saveas(fig,['D:\Sub_beta_up',num2str(sss),'_c_',num2str(cl),'_t_',num2str(tim),'_f_',num2str(fre),'__'],'png')
                print('-depsc2', ['I:\Sub_',num2str(sss),'_c_',num2str(cl),'_t_',num2str(tim),'_f_',num2str(fre),'__.eps'])
                %                 view([-90 90])
                close
            end
        end
    end
end