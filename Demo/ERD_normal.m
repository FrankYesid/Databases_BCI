
load('G:\Dropbox\ERD\Codes\TP\Matlab_con_csp\BCICIV_1.mat')

%%
X_suj = cell(7,1);
fs = 100;
for s = [3,8,9]
    ns = numel(X);
    N_channels = size(X{1}{1},2);
    
    % %%%%%%%%%%%%%%%%%%%%%%%%% ERDs spectrogram %%%%%%%%%%%%%%%%%%%%%%%%% %% % parametros de la STFT
    % Nw = 1;                                       % tama?o de la ventana en seg
    Nw = 0.5*fs;                                    % tamaño de la ventana.
    window = hamming(Nw);                           % Vector de la ventana.
    noverlap = floor(0.90*Nw);                      % interpolación de la ventana.
    X_suj{s}=cell(2,1);
    % celdas almacenadoras de spectrogram.
    for Clase = 1:2                                  % clase o clases seleccionadas.
        su = s ;                                   % calculo del espectrograma relacionada a EEG.
        Sujeto = X{su};                                                         % sujeto seleccionado.
        etiqueta = y{su};
        % etiqueta seleccionada.
        % selecciona los trail de las etiquetas solicitadas.
        % Izq_Der = Sujeto(ismember(etiqueta,[1 2]));
        Izq_Der = Sujeto(ismember(etiqueta,Clase));                             % seleccion de los trails segun la clase o clases seleccionadas.
        % escojer el que tenga la menor cantidad de trials
        N_trial = length(Izq_Der);                                              % cantidad de trials.
        Canales = [1:N_channels]; N_canal = length(Canales);                            % cantidad de canales.
        % Inicializar matrices
        temp = spectrogram(Izq_Der{1}(:,1), window, noverlap,Nw);               % tamaño temporal de los datos en cada sujeto.
        temp = size(temp,2);
        X_suj{su}(Clase) = {zeros(N_trial,N_canal,floor((Nw/2)+1),temp)};              % Trials x canal x frecuencia x tiempo.
        for tri = 1:N_trial
            fprintf(['sujeto: ' num2str(su) ' ...trial: ' num2str(tri) ' de ' num2str(N_trial) '\n'])
            for cnl = 1:N_canal
                % signal
                signal = Izq_Der{tri}(:,Canales(cnl));
                
                % Calcular STFT
                [X_Class, f, t] = spectrogram(signal, window, noverlap,Nw,fs);
                X_Class = abs(X_Class);                                         % absoluto para tener la señal (freq,tiempo)
                
                % Almacenar espectrogramas
                X_suj{su}{Clase}(tri,cnl,:,:) = X_Class;                               % Trials x canal x frecuencia x tiempo
            end
        end 
    end
    % tiempo de referencia.
    %     save(['X_suj_' num2str(s) '.mat'],'X_suj')
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% ERD - Normal %%%%%%%%%%%%%%%%%%%%%%%%%%% %% % calculo del ERD
t1 = 0.5;                                         % tiempo inicial de referencia.
t2 = 2;                                         % tiempo final de referencia.
temp1 = abs(t - t1);
min1 = min(temp1);
temp2 = abs(t - t2);
min2 = min(temp2);
ul = find(temp1 == min1);
up = find(temp2 == min2);
% r_nc = cell(ns,1);                               % celda para esperanza respecto al tiempo de referencia.
% r_c = cell(ns,1);                                % celda para esperanza respecto a los trials (intentos).
% m_c = cell(ns,1);                                % celda para esperanza respecto a los trials (intentos).
ERD = cell(ns,1);                                % celda para el ERD.

%% mean ERD
for s = [1,2,6,7]
    %     load(['X_suj_' num2str(s) '.mat'])
    fprintf(['Subject ' num2str(s) '\n'])
    for cl =1:2
        r_nc = squeeze(mean(X_suj{s}{cl}(:,:,:,ul:up),4));
        r_c = squeeze(mean(r_nc,1));
        m_c = squeeze(mean(X_suj{s}{cl},1));
        ERD{s}{cl} = bsxfun(@times,m_c,1./r_c) - 1;
    end
end



%% %%%%%%%%%%%% Headplot for channel ciclo spectrogram - 3D %%%%%%%%%%% %%%
set(0,'DefaultFigureWindowStyle','docked')
load('electrodesBCICIV1.mat')
% load('ERD_normal.mat')
% grafica cada figura en una misma pantalla.
% posición de las graficas
pos_ = [4,6,11,12,13,14,15,16,17,20,21,22,23,24,25,26,28,29,30,31,33,34,...
    35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,51,52,53,54,56,57,58,...
    59,60,61,62,65,66,67,68,69,70,71,76,78,85,87];  % ubicación de cada una de las graficas.
% Rango de freq.
% frang = [5 16];                                                             % rango de frecuencias a graficar.
% Rango de tiempo
% trang = [5 595];                                                           % rango de tiempo a graficar.
% rango de graficar limites
rt = [1 5]; % tiempo                                                        % rango en el plot del tiempo.
rf = [5 30];% freq                                                          % rango en el plot de las frecuencias.
% graficador en 3D todos los canales
for sub = [2,7]
    for cl = 1:2
        figure
        for ch = 1:59
            subplot(10,9,pos_(ch))
            %             surf(t,f,squeeze(ERD{sub}{cl}(ch,:,:)))
            imagesc(t,f,squeeze(ERD{sub}{cl}(ch,:,:)))
            title(electrodes(ch))
            shading interp
            %             view(-18,38)
            xlim([0 8])
            %             view(0,90)
            yticks([])
            xticks([])
            colormap jet
            grid on
            hold on
            line([2,2],[1,50],'LineWidth',2,'color','r')
            line([4,4],[1,50],'LineWidth',2,'color','r')
            hold off
            drawnow
        end
        colorbar
        suptitle(['Subject ' num2str(sub) ' Class ' num2str(cl)])
    end
end