%% main topolot
% Limpiar la memoria almacenada, graficas creadas y la ventana de comandos.
clear all; close all; clc

% Lugar de almacenamiento de los topoplots.
save = '';

% Lugar de funciones adicionales si son necesarias.

%  1- Hacer un topoplot general que te permita identificar todos los
%  canales de 128 posiciones, y permita la interacción de diferentes
%  herramientas como archivos .mat, .edf, .loc.

% Lectura de los canales que se quieren usar.
load();
sig = labels; % coloca la variable que contenga los nombres de los canales LABELS.

% Cargar la configuración de los canales del topoplot.


% Cargar el tipo de la cabeza.
% 2D
load('HeadModel.mat')   % model of the head.

% Cargar la relevancia de los diferentes canales.
rel = [];
% Si es necesario coloca las relevancias entre [0-1];

% Sistema comparador de los canales usados por los labels utilizados
%  = 



% Función encargada de la creación del topoplot
MyTopo_fun(rel,sel,M1.xy,M1.lab,[min(rel) max(rel)],0,0,t_e,database,HeadModel,type,tex,c)
% caxis([-1,1])
% colormap(gca,cmap);
title(['Sujeto ' num2str(s)],'Interpreter','latex')
axis square

