load('BCICIV_2a.mat');                          % carga base de datos EEG.
name = 'BCICIV_2a.mat';
fprintf('Loaded data: %s\n',name);

ns = numel(X);

[~,N_channel] = size(X{1}{1});
sub = 2;
Xd_ = cell(1,ns);
lb = cell(1,ns);
labels = [1 2];
N_class = 2;
for i = 1:ns
    ind{i} = ismember(y{i},labels);
    Xd_{i} = X{i}(ind{i});                      %belonging to labels
    lb{i} = y{i}(ind{i});
end

sig = cell(1,ns);
for s = sub
    for class = 1:N_class
        Xr_ =  Xd_{s}(ismember(lb{s},class));
        for  tr = 1:numel(Xr_)
            signal = Xr_{tr};
            sig_cov = sqrt(trace(cov(signal)));
            data(tr,:,:) = Xr_{tr}./sig_cov;
%                         data(tr,:,:) = signal;
        end
%                 sig{s}{class}=sig_cov;
        sig{s}{class}=data;
    end
end

Nw = fs;                                        % tamaño de la ventana.
window = hamming(Nw);                           % Vector de la ventana.
noverlap = floor(0.98*Nw);                      % interpolación de la ventana.
X_suj = cell(N_class,1);                              % celdas almacenadoras de spectrogram.

for s = sub
    for Clase = 1:N_class                           % clase o clases seleccionadas.                             % etiqueta seleccionada.
        X_se = sig{s}{Clase};                                     % seleccion de los trails segun la clase o clases seleccionadas.
        N_trial = size(X_se,1);                     % cantidad de trials.
        for tri = 1:N_trial
            fprintf(['sujeto: ' num2str(s) ' de 9' ' Class: ' num2str(Clase) ' ...trial: ' num2str(tri) ' de ' num2str(N_trial) '\n'])
            signaltr = squeeze(X_se(tri,:,:));
            for cnl = 1:N_channel
                signal = signaltr(:,cnl);
                [X_Class, fre, time] = spectrogram(signal, window, noverlap,Nw,fs);% Calcular STFT
                X_Class = abs(X_Class);             % absoluto para tener la señal (freq,tiempo)
                d(cnl,:,:) = X_Class;               % Trials x canal x frecuencia x tiempo % Almacenar espectrogramas.
            end
            dat(tri,:,:,:) = d;
        end
        X_suj{s}{Clase} = dat;
    end
end

%% %%%%%%%%%%%% Headplot for channel ciclo spectrogram - 2D %%%%%%%%%%% %%
set(0,'DefaultFigureWindowStyle','docked')                                  % grafica cada figura en una misma pantalla.
% posición de las graficas
posi = [4 9 10 11 12 13 15 16 17 18 19 20 21 23 24 25 26 27 31 32 33 39];   % ubicación de cada una de las graficas.

for s = sub
    i = 1;
    for class = 1:N_class
        figu = figure;
        for ch = 1:22
            subplot(6,7,posi(ch))
            surf(time,fre,squeeze(std(X_suj{s}{class}(:,ch,:,:))))
            %             histogram(sig{s}{class},'BinLimits',[20 100])
            %             set(gca,'ColorScale','log')
            title(['Ch: ' num2str(ch)])
            shading interp
            %             set(gca,'ZScale','log')
            %             zlim([0 1])
            %             set(gca,'YScale','log')
            %                         caxis([0 1])
            view([0 90])
            colorbar
            if ch == 22
                xlabel('time (seg)')
            end
            %         ylim([1 40])
            xlim([1 time(end)])
            ylim([1 30]);
            colormap jet
            %             xlim([20 100])
            drawnow
            i = i+1;
        end
        %         suptitle(['Sujeto ' num2str(sub) ' clase ' num2str(clas)])
        suptitle(['Sujeto ' num2str(s) ' - clase ' num2str(class)])
        %         saveas(figu,['C:\Users\frany\Desktop\Frank\S' num2str(sub) 'c' num2str(class) 'std_.fig']);
        saveas(figu,['C:\Users\Lufe\Desktop\Dispersion graficas\s' num2str(s) 'c' num2str(class) 'std.fig']);
        close all
    end
end
