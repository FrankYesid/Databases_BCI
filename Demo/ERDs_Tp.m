% clear; clc
% load('G:\Dropbox\ERD\Codes\TP\Matlab_con_csp\BCICIV_1.mat')
%load('F:\Database\BCICIV_2a\BCICIV_2a.mat')
% load('BCIIV1_TPbyclass.mat')
%
% % ERDs
% for s = [2,7]
%     for class = [1,2]
%         dats = X{s}(y{s}==class);
%
%         for fr = 1:numel(Wt{s})
%             tic
%             W_tp = Wt{s}{fr}{class};
%             for tr = 1:numel(dats)
%                 signals = dats{tr};
%                 N_channel = size(signals,2);
%                 parfor ch = 1:N_channel
%                     E(:,ch) = conv(signals(:,ch),W_tp(ch,:),'same');
%                 end
%                 ERD_t{s}{class}(tr,fr,:,:) = E;
%             end
%             fprintf(['ERD-TP: Subject ' num2str(s) ' Freq: ' num2str(fr) ' Trial: ' num2str(tr) ' time: ' num2str(toc) '\n'])
%         end
%         ERDs{s}{class} = squeeze(mean(ERD_t{s}{class},1));
%     end
% end
% save('ERDs_tp.mat','ERDs')
% %% plot
% load('ERDs_tp.mat','ERDs')
% load('electrodesBCICIV1.mat')
% %
% pos_ = [4,6,11,12,13,14,15,16,17,20,21,22,23,24,25,26,28,29,30,31,33,34,...
%     35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,51,52,53,54,56,57,58,...
%     59,60,61,62,65,66,67,68,69,70,71,76,78,85,87];
% set(0,'DefaultFigureWindowStyle','docked')
% for s = [2,7]
%     for class = 1:2
%         fig=figure;
%         set(fig,'name',['Subject ' num2str(s)])
%         for ch = 1:59
%             subplot(10,9,pos_(ch))
%             imagesc(0:1/100:8-1/100,1:39,ERDs{s}{class}...
%                 (:,:,ch),[-10 10])
%             grid on
%             yticks([])
%             xticks([])
%             colormap jet
%             hold on
%             line([2,2],[1,39],'LineWidth',2,'color','r')
%             line([4,4],[1,39],'LineWidth',2,'color','r')
%             hold off
%             title(electrodes(ch))
%             xlim([0 7])
%             if ch >= 58
%                 yticks([1 39])
%                 xticks([0 2 4 7])
%             end
%             drawnow
%         end
%         colorbar
%         suptitle(['Subject ' num2str(s) ' Class ' num2str(class)])
%     end
% end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear; clc
% load('G:\Dropbox\ERD\Codes\TP\Matlab_con_csp\BCICIV_1.mat')
% load('BCIIV1_TPbyclass4-30.mat')
load('F:\Database\BCICIV_2a\BCICIV_2a.mat')
% X_suj = cell(9,1);
% fs = 100;
X_suj = cell(numel(X),1);
% ERDs
clc

for s = [8]
    ns = numel(X);
    N_channels = size(X{s}{1},2);
    
    % %%%%%%%%%%%%%%%%%%%%%%%%% ERDs spectrogram %%%%%%%%%%%%%%%%%%%%%%%%% %% % parametros de la STFT
    % Nw = 1;                                       % tama?o de la ventana en seg
    Nw = 0.5*fs;                                    % tamaño de la ventana.
    window = hamming(Nw);                           % Vector de la ventana.
    noverlap = floor(0.90*Nw);                      % interpolación de la ventana.
    X_suj{s}=cell(2,1);
    
    for Clase = 1:2
        dats = X{s}(y{s}==Clase);
        tic
        %         W_tp = Wt{s}{1}{Clase};
        su = s;                                    % calculo del espectrograma relacionada a EEG.
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
        %         X_suj{su}(Clase) = {zeros(N_trial,N_canal,floor((Nw/2)+1),temp)};             % Trials x canal x frecuencia x tiempo.
        E=zeros(N_channels,floor((Nw/2)+1),temp);
        for tr = 1:numel(dats)
            signals = dats{tr};
            N_channel = size(signals,2);
            tic
            for ch = 1:N_channel
                if Ch_img{s}{1}{Clase}(1,ch) == 1
                    continue
                else
                    %                     ttem = conv(signals(:,ch),W_tp(ch,:),'same');
                    [X_Class, f, t] = spectrogram(signals(:,ch), window, noverlap,Nw,fs);
                    X_Class = abs(X_Class);
                    ta = size(X_Class,1);
                    if ta <= 26
                        E(ch,1:26,:) = X_Class;
                    else
                        E(ch,:,:) = X_Class;
                    end
                end
            end
            ERD_t{s}{Clase}(tr,:,:,:) = E;
            fprintf(['ERD-TP: Subject ' num2str(s) ' Trial: ' num2str(tr) ' time: ' num2str(toc) '\n'])
        end
        
        ERDs{s}{Clase} = squeeze(mean(ERD_t{s}{Clase},1));
    end
end

% save('ERDs_tp_2.mat','ERDs')
%% plot
% load('ERDs_tp_2.mat','ERDs')
load('electrodesBCICIV1.mat')
%
pos_ = [4,6,11,12,13,14,15,16,17,20,21,22,23,24,25,26,28,29,30,31,33,34,...
    35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,51,52,53,54,56,57,58,...
    59,60,61,62,65,66,67,68,69,70,71,76,78,85,87];
set(0,'DefaultFigureWindowStyle','docked')
for s = [2,7]
    for Clase = 1:2
        fig=figure;
        set(fig,'name',['Subject ' num2str(s)])
        for ch = 1:59
            if Ch_img{s}{1}{Clase}(ch) == 1
                continue
            else
                subplot(10,9,pos_(ch))
                imagesc(t,f,squeeze(ERDs{s}{Clase}(ch,:,:)))
                grid on
                %             view(0,90)
                %                          yticks([0 50])
                %             xticks([0 2 4 7])
                yticks([])
                xticks([])
                colormap jet
                hold on
                line([2,2],[1,59],'LineWidth',2,'color','r')
                line([4,4],[1,59],'LineWidth',2,'color','r')
                hold off
                title(electrodes(ch))
                xlim([0 7])
                %             ylim([0 15])
                %             if ch >= 58
                %                 yticks([1 39])
                %                 xticks([0 2 4 7])
                %             end
                drawnow
                grid on
            end
        end
        colorbar
        suptitle(['Subject ' num2str(s) ' Class ' num2str(Clase)])
    end
end
