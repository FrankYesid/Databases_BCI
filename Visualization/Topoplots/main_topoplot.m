%% MAIN TOPOPLOT - SELECTION OF TOPOPLOT
% load('mycamp')    % selección determinada de colores.
% clear all; close all; clc

% Direccion de la base de datos
% SUBJECTS_DIR = 'D:\Luisa\Dropbox\ERD\results_ERDfc_subjects\BCI'; % Ubicacion de la Información de cada sujeto.
SUBJECTS_DIR = 'D:\Dropbox\ERD\results_ERDfc_subjects\BCI'; % Ubicacion de la Información de cada sujeto.
% Seleccion de la base de datos de BCI.
% (1) BCICIV_1, (2) BCICIV_2a_, (3) GIGASCIENCE_, (4) BCIIII_4a_
database = 2;
COHORT   = 'BCICIV_2a_';
SUBJECTS = dir([SUBJECTS_DIR filesep '*' COHORT '*']);
SUBJECTS = struct2cell(SUBJECTS);
SUBJECTS = SUBJECTS(1,:)';

% Sujetos Seleccionados
subjects = [5];
% Seleccion del tipo de topoplot requerido.
% (1) Selection of channels.
% (2) Selection of channels with relevance.
% (3) Plots of signals en each electrode.
% (4) Creating data mask.
type   = 3;
% Las figuras en subplots (0) OFF (1) ON
subpl  = 0;
% Para colocar los nombres de cada uno de los canales.
tex    = 1;

% Colocar imagenes de una misma ventana.
set(0,'DefaultFigureWindowStyle','docked')

% fold seleccionado para el W.
fold   = 5;
% Caracterististica seleccionada.
caract = 1;

% Colores de la grafica
color1 = [014,041,075]/255;         % azul oscuro
color2 = [1 1 1];                   % blanco
% color3 = [0.9763,0.9831,0.0538];  % Naranja claro brillante
color3 = [236 124 038]/255;         % Naranja intenso
% Construye el colormap
c = [[linspace(color1(1),color2(1),32),linspace(color2(1),color3(1),32)]',...
    [linspace(color1(2),color2(2),32),linspace(color2(2),color3(2),32)]',...
    [linspace(color1(3),color2(3),32),linspace(color2(3),color3(3),32)]'];

% Carga los pesos del laso.
% load('G:\Dropbox\ERD\results_ERDfc_subjects\BCI\BCICIV_2a\pesos_lasso_prom_frec.mat')
% cargar vector de relevancia de canales.
% load('')

% definir parametros de filter bank
f_low  = 4;
f_high = 40;
Window = 4;
Ovrlap = 2;
filter_bank = [f_low:Ovrlap:f_high-Window;...
    f_low+Window:Ovrlap:f_high]';
orden_filter = 5;
labels = [1 2];

for s = subjects
    figure(s)
    % posicion inicial del subplot.
    a      = 1;
    %     load([SUBJECTS_DIR filesep SUBJECTS{s} filesep 'results' filesep 'ERDCSP_rho.mat']) % Carga del rho y del W.
    %     load([SUBJECTS_DIR filesep SUBJECTS{s} filesep 'results' filesep 'Ex_Giga_2W_vDic.mat'])
    rho = zeros(1,22);
    if sum(isnan(rho))>1; rho(isnan(rho)==1)=1;  end                                 % para los nan detectados
    if s ==2; bb = [1,2,4,7,16]; end
    if s ==5; bb = [1,2,5,8,11:14]; end
    if s ==8; bb = [1:5,8,13]; end
    for band = bb                       % bandas seleccionadas.
        switch database                 % database seleccionada.
            
            case 1
                load('BCICIV_1\electrodesBCICIV1.mat') % se encuentra la ubicacion de los electrodos y sus nombres.
                load('HeadModel.mat')   % model of the head.
                Fil = 2; Col = 2;       % topoplot en subplots para cada uno de los sujetos.
                t_e = 70;               % tamaño visual de los electrodos seleccionados.
                sel = (1:59);           % Orden de los canales segun una relevancia.
                M1.xy = pos(:,1:2);     % posicion de los canales.
                M1.lab = electrodes;    % nombre de los canales.
                rel = linspace(1,10,numel(M1.lab));     % para colocar mas importancia al electrodo.
                %                 sub = [59,59,59,59];    % cantidad de canales seleccionados para cada sujeto.
                
            case 2
                load('BCICIV_2a\electrodesBCICIV2a.mat')
                %load('BCICIV_2a\labels.mat')
                %load('BCICIV_2a\layout.mat')
                load('HeadModel.mat')   % model of the head.
                Fil = 3; Col = 3;             % topoplot en subplots para cada uno de los sujetos.
                t_e = 70;               % tamaño visual de los electrodos seleccionados.
                M1.xy = elec_pos(:,1:2);% posicion de los canales.
                M1.lab = Channels;      % nombre de los canales.
                mask = reshape(rho,[17 22]) <= thresh; %
                sel = find(mask(band,:)==1); % Orden de los canales segun una relevancia.
                % relevancia para los canales la variable es w.
                tmp = W{fold,band}(caract,:); % rel = linspace(1,10,numel(M1.lab));     % para colocar mas importancia al electrodo.
                rel = zeros(1,22);
                rel(sel) = tmp;
                
                % sub = [22,22,22,22];    % cantidad de canales seleccionados para cada sujeto.
                
            case 3
                load('Gigasc\posiciones.mat')
                load('HeadModel.mat')   % model of the head.
                t_e = 70;               % tamaño visual de los electrodos seleccionados.
                sel = (1:64);           % Orden de los canales segun una relevancia.
                M1.xy = floor(M1.xy(:,[1 2]).*100); % posicion de los canales.
                M1.xy(28,1) = -90; % posicion de los canales.
                M1.lab = M1.lab;        % nombre de los canales.
                rel = linspace(1,10,numel(M1.lab));    % para colocar mas importancia al electrodo.
                %                 sub = [64,64,64,64];    % cantidad de canales seleccionados para cada sujeto.
                
            case 4
                load('BCICIIIIVa\electrodesBCICIIIIVa.mat')
                load('HeadModel.mat')   % model of the head.
                t_e = 100;              % tamaño visual de los electrodos seleccionados.
                sel = (1:118);          % Orden de los canales segun una relevancia.
                M1.xy = pos;            % posicion de los canales.
                M1.lab = electrodes;    % nombre de los canales.
                rel = linspace(1,10,numel(M1.lab));    % para colocar mas importancia al electrodo.
                %                 sub = [118,118,118,118];% cantidad de canales seleccionados para cada sujeto.
                
            case 5
                
        end
        % end
        
        % for s = subjects
        if subpl == 1
            rel = 1:22;
            subplot(1,numel(sub),a)
            MyTopo_fun(rel,sel(1:sub(a)),M1.xy,M1.lab,[min(rel) max(rel)],0,0,t_e,database,HeadModel,type,tex,c)
            % caxis([-1,1])
            % colormap(gca,cmap);
            title(['Sujeto ' num2str(s)],'Interpreter','latex')
            axis square
            axis off
            a=a+1;
        else
            %             subplot(2,4,a)
            rel = rel/max(abs(rel));
            MyTopo_fun(rel,sel,M1.xy,M1.lab,[min(rel) max(rel)],0,0,t_e,database,HeadModel,type,tex,c)
            % caxis([-1,1])
            % colormap(gca,cmap);
            title(['Sujeto ' num2str(s)],'Interpreter','latex')
            axis square
            
            axis off
            title([num2str(tmpuF(s,band)) ' - Banda: ' num2str(filter_bank(band,:))])
            a=a+1;
        end
    end
    %     fig = gca;
    %     saveas(fig,['G:\Dropbox\ERD\results_ERDfc_subjects\BCI\BCICIV_2a' filesep 'Sub' num2str(s) 'Topoplot.fig'])
    %     saveas(fig,['G:\Dropbox\ERD\results_ERDfc_subjects\BCI\BCICIV_2a' filesep 'Sub' num2str(s) 'Topoplot.png'])
    %     saveas(fig,['G:\Dropbox\ERD\results_ERDfc_subjects\BCI\BCICIV_2a' filesep 'Sub' num2str(s) 'Topoplot'],'epsc')
end

