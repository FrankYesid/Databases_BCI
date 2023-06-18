%% MAIN TOPOPLOT - CONNECTIVITY TOPOPLOT
clear; close all; clc

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
subjects = 5;
% Seleccion del tipo de topoplot requerido.
% (1) Selection of channels.
% (2) Selection of channels with relevance.
% (3) Plots of signals en each electrode.
% (4) Creating data mask.
type   = 4;
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

load('BCICIV_2a\electrodesBCICIV2a.mat')
% load('BCICIV_2a\labels.mat')
% load('BCICIV_2a\layout.mat')
load('HeadModel.mat')    % model of the head.
Fil = 3; Col = 3;              % topoplot en subplots para cada uno de los sujetos.
t_e = 20;                        % tamaño visual de los electrodos seleccionados.
M1.xy = elec_pos(:,1:2);% posicion de los canales.
M1.lab = Channels;        % nombre de los canales.
sel = 1:22;                     % canales seleccionados.
rel = ones(1,22);            % relevancia de los canales.

for s = subjects
    figure(s)
    pos = M1.xy;
    label = M1.lab;
    
    % database = 3;
    % tam - tamaño de la cabeza y posición de los electrodos.
    temm = 0.1;
    r_tam = 0.5;
    tam = 1; tam2 = 0.9-temm; tam3 = 0.99-temm;
    for i = 1:2
        pos(:,i) = tam2.*((pos(:,i)-min(pos(:,i)))/(range(pos(:,i)))-0.5);
    end
    xc = HeadModel(1,:);
    yc = HeadModel(2,:);
    hold on
    % Topoplot
    x = pos(:,1)+pos(1,1);
    y = pos(:,2)+pos(1,2);
    tmp = [x,y,x*0]*rotz(2);
    x = tmp(:,1); y = tmp(:,2);
    pos(:,1) = x; pos(:,2) = y;
    GS = 800;
    xi = linspace(min(x)-0.9, max(x)+0.9, GS);       % x-axis for interpolation (row vector)
    yi = linspace(min(y)-0.9, max(y)+0.9, GS);       % y-axis for interpolation (row vector)
    Method_grid ='v4'; %'natural';
    % 'linear'  Triangle-based linear interpolation (default)
    % 'cubic' Triangle-based cubic interpolation
    % 'nearest' Nearest neighbor interpolation
    % 'v4' MATLAB 4 griddata method
    [Xi, Yi, Zi] = griddata([x' xc], [y' yc], [rel(:);zeros(numel(xc),1)], xi', yi,Method_grid); % interpolate the topographic data
    %% Creating data mask
    [TH,R] = cart2pol(Xi,Yi);
    Zi(R>r_tam) = NaN;
    deltax = xi(2)-xi(1); % length of grid entry
    deltay = yi(2)-yi(1); % length of grid entry
    h = surf(Xi-deltax/2, Yi-deltay/2+0.002, zeros(size(Zi)), Zi,'EdgeColor', 'none', 'FaceColor', 'flat');hold on
    %     shading interp % interpola colores.
    scatter(tam3.*x(sel),(tam3+0.04).*y(sel),t_e,'b','filled')
    %     scatter(tam*xc(sel),tam*yc(sel)+0.001,t_e,'b','filled')
    % text(pos(sel,1)*tam3-0.02,pos(sel,2)*tam3+0.02,label(sel),'Interpreter',...
    %         'latex','ButtonDownFcn',{@lineCallback,database},'Color','black',...
    %         'FontSize',10);
    %     plot(tam*xc+pos(10,1),tam*yc+0.001+pos(10,2),'k','LineWidth',3)
    plot(tam*xc,tam*yc+0.001,'k','LineWidth',3)
    caxis([0 1])
    axis off
    %         hold off
    axis square
    %     view([-90 90])
    title(['Sujeto ' num2str(s)],'Interpreter','latex')
end


