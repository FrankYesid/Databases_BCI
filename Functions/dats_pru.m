clear all; close all; clc
%%

% % Direccion de la base de datos
% SUBJECTS_DIR = 'D:\BCI';
SUBJECTS_DIR = 'F:\BCI';
% % Direccion del fold de las funciones
% % addpath(genpath('C:\Users\lfvelasquezm\Dropbox\ERD\Codes\TP\Matlab_wang\csp\CSP_fun\functions'))
% addpath(genpath('C:\Users\lfvelasquezm\Desktop\frank\functions'))
% addpath(genpath('C:\Users\frany\Dropbox\Event-related\Codes\TP\Matlab_wang\csp\CSP_fun\functions'));
addpath(genpath('G:\Dropbox\ERD\Codes\TP\Matlab_wang\csp\CSP_fun\functions'))
%
% %% DataBase
% % BCIIII_4a_
% % BCICIV_2a_
% % GIGASCIENCE_
%
COHORT = 'BCICIV_2a_';
SUBJECTS = dir([SUBJECTS_DIR filesep '*' COHORT '*']);
SUBJECTS = struct2cell(SUBJECTS);
SUBJECTS = SUBJECTS(1,:)';
%
%
% %% grilla de busqueda
param = linspace(0,0.9,100);

experiment_name = 'Ex_Giga_vLFVM';
experiment_name1 = 'Ex_Giga_vLFVM_mod';
% SS = [37 32 12 18 42 34 3 7 35 33 21 2 4 39 29 43 28]; % UNO BUENO Y UNO MALO%%%%INDEXACDOS DE ACIERDO A CSP
% SS = [37 32 12 18 42 34 3 7];
SS = 1:9;
% if strcmp(COHORT,'GIGASCIENCE_')
%     SubInd = [50,14];
%     SS(SubInd) = [];
% end

%% paramaters definition
Rep = 10000;  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% numero de repeticiones
PPval = cell(numel(SS),1);
rho = [];
tstart = 0;
tend = 2;
% definir parametros de filter bank
f_low  = 4;
f_high = 40;
Window = 4;
Ovrlap = 2;
filter_bank = [f_low:Ovrlap:f_high-Window;...
    f_low+Window:Ovrlap:f_high]';
orden_filter = 5;
labels = [1 2];

for s = SS
    load(['G:\Dropbox\ERD\results_ERDfc_subjects\BCI\' filesep SUBJECTS{s} filesep experiment_name 'rho.mat'],'rho_s')
    if sum(isnan(rho_s))>1; rho_s(isnan(rho_s)==1)=1;  end %para los nan detectados
    threshold = linspace(min(rho_s),max(rho_s),100);    
    load(['G:\Dropbox\ERD\results_ERDfc_subjects\BCI\' filesep SUBJECTS{s} filesep experiment_name1 'acc_mod_lu.mat'],'acc','table');
    macc = nanmean(acc,1);
    mmacc(s) = max(macc(:));
    table1(s,:) = param_subject(acc,threshold,param);
end