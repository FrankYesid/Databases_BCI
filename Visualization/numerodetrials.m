clear all; close all; clc
%%

% % Direccion de la base de datos
% SUBJECTS_DIR = 'D:\BCI';
SUBJECTS_DIR = 'F:\BCI';
% % Direccion del fold de las funciones
% % addpath(genpath('C:\Users\lfvelasquezm\Dropbox\ERD\Codes\TP\Matlab_wang\csp\CSP_fun\functions'))
addpath(genpath('C:\Users\lfvelasquezm\Desktop\frank\functions'))
% addpath(genpath('C:\Users\frany\Dropbox\Event-related\Codes\TP\Matlab_wang\csp\CSP_fun\functions'));
%
% %% DataBase
% % BCIIII_4a_
% % BCICIV_2a_
% % GIGASCIENCE_
%
COHORT = 'GIGASCIENCE_';
SUBJECTS = dir([SUBJECTS_DIR filesep '*' COHORT '*']);
SUBJECTS = struct2cell(SUBJECTS);
SUBJECTS = SUBJECTS(1,:)';
%
%
% %% grilla de busqueda
param = linspace(0,0.9,100);

experiment_name = mfilename;
%
% SS = [37 32 12 18 42 34 3 7 35 33 21 2 4 39 29 43 28]; % UNO BUENO Y UNO MALO%%%%INDEXACDOS DE ACIERDO A CSP
SS = 1:52;
% if strcmp(COHORT,'GIGASCIENCE_')
%     SubInd = [7,9,29,32,46,49,26,31];
%     SS(SubInd) = [];
% end
labels = [1 2];

%
for s = SS
    %     clearvars -except s SS rho experiment_name COHORT param SUBJECTS SUBJECTS_DIR Acc table PPval Rep tstart tend
    %
    load([SUBJECTS_DIR filesep SUBJECTS{s} filesep 'eeg' filesep 'raw.mat'])
    y = y(:);
    ind = ismember(y,labels);
    y = y(ind);
    X = X(ind);
    
    trials(s) = numel(X); 
    fprintf(['Sujeto: ' num2str(s) '...' num2str(trials(s)) '\n'])
end
