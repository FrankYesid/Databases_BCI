%% Database
clear all; clc
%% load data
% load('BCICIV_2a.mat');
load('D:\MATLAB\CKA_BCI\Old\BCICIV_2a.mat') % load files
%% ERD con Wavelet band

t = 0:1/250:7-(1/250);
t1 = 0;
t2 = 2;
lv = 5;  %  filter order
% sujeto = 3;
% name = 'ad3';
for su = 9
    %     su = %1:9
    
    Sujeto = X{su};
    etiqueta = y{su};
    
    % escojer los train de las etiquetas 1 y 2
    Data = Sujeto(ismember(etiqueta,[1 2]));
    
    % escojer el que tenga la menor cantidad de trials
    N_trial = length(Data);
    Canales = [1:22]; N_canal = length(Canales);
    
    for tri = 1:N_trial
        fprintf(['sujeto: ' num2str(su) ' de ' '9' ' ...trial: ' num2str(tri) ' de ' num2str(N_trial) '\n'])
        % signal
        signal = Data{tri}(:,:);
        % function filterwavelet
        [X_dat, T] = FilterWave(signal,'sym2',lv);
        Y_suj{su}(tri,:,:,:) = X_dat;
    end
    
    
    %%
    % t1 = 0;
    % t2 = 2;
    fprintf(['ERD = m_c/r_c\n'])
    % for
    fprintf(['sujeto: ' num2str(su) ' de ' '9' '\n'])
    % ERD = m_c/r_c
    ERD= ERDS(Y_suj{su},t1,t2);
    % end
%     erders{su} = ERD;
%     save ERER ERD
end

load('ERD.mat')
erders{9}= ERD;
save ERD erders


