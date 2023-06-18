%% Programa Wang2012
%-------------------------------------------------------------------------
% Inputs:
% X{subject}{trial}(samples x chan): Data to analysis
% y{subject}{trial}
% fs: Sample frequency
% labels{channels}
%-------------------------------------------------------------------------
% Outputs:
% Xf: Filtered data
% Xf{trial}(samples x chan)
% segun el uso de spectrogram o wavelet packet (con maximo o con entropia).
% ERD{subject}(chan x freq x time)
% SampEnc1 ans SampEnc2 -> Sample entropy
% ttest(SampEnc1,SampEnc2) with Default is 0.05 for 5% significance.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2018 Signal Processing and Recognition Group
% F.Y. Zapata C.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%% Seccion A -wang %%%%%%%%%%%%%%%%%%%%%%%%%% %%
% clear all; clc                                  % limpia la ventada de comandos y también la memoria del espacio de trabajo.
% %% Manejo del codigo en general
% database = 2;                                   % Database 1-> BCICIV_2a or 3-> BCIIII_4a.
% norm = 0;                                       % Normalization 1->ON 0->off.
t1 = 0;                                         % tiempo inicial de referencia.
t2 = 2;                                         % tiempo final de referencia.
t3 = 2.5;                                       % tiempo inicial de ventana deslizante
t4 = 4.5;                                       % tiempo final de ventana deslizante.
spec = 0;                                       % spectrogram de la señal en cada uno de los canales
filt = 1;                                       % filter [2 40]
% % sign_lv = 0.05;                                 % significance level.
% %% Load Database
% if database == 1                                % carga base de datos EEG. BCICIV_2a.mat
%     load('F:\Database\BCICIV_2a.mat');
%     name = 'BCICIV_2a';
% elseif database == 2
%     load('F:\Database\BCICIV_1\dataBase.mat')   % carga base de datos EEG. BCICIV_1
%     fs = 100;
%     name = 'BCIIII_4a';
% elseif database == 3                            % carga base de datos EEG. BCIIII_4a.mat
%     X_ = cell(5,1);
%     y_ = cell(5,1);
%     for s = 1:5
%         load(['F:\BCI\BCIIII_4a_0' num2str(s) '\eeg\raw.mat'])
%         X_{s,1} = X;
%         y_{s,1} = y;
%         fs = fs;
%     end
%     X = X_; clear X_
%     y = y_; clear y_
%     name = 'BCIIII_4a';
% end
% % load('labels.mat');                           % carga nombres de los canales.
% % load('layout');                               % carga localizacion.
% fprintf('Loaded data: %s\n',name);
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%% Normalization %%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% if norm == 1
%     X = cellfun(@(Y) cellfun(@(X) (X - mean(X))/ std(X),Y,...
%         'UniformOutput',false),X,'UniformOutput',false);
% end
% %% %%%%%%%%%%%%%%%%%%%%%%%%% ERDs-spectrogram %%%%%%%%%%%%%%%%%%%%%%%%% %% % parametros de la STFT
% if spec == 1
ns = numel(X);                                  % numero de sujetos.
% Nw = 1;                                       % tamaño de la ventana en seg
Nw = 0.8*fs;                                    % tamaño de la ventana.
window = hamming(Nw);                           % Vector de la ventana.
noverlap = floor(0.95*Nw);                      % interpolación de la ventana.
n_class = 2;
X_suj = cell(n_class,1);
for clase = [1 2]                               % clase o clases seleccionadas.
    X_suj{clase} = cell(ns,1);                  % celdas almacenadoras de spectrogram.
    fprintf(['Clase: ' num2str(clase) '\n'])
    for su = 1:ns                               % calculo del espectrograma relacionada a EEG.
        tic
        Sujeto = X{su};                         % sujeto seleccionado.
        etiqueta = y{su};                       % etiqueta seleccionada.
        % Izq_Der = Sujeto(ismember(etiqueta,[1 2]));                           % selecciona los trail de las etiquetas solicitadas.
        Izq_Der = Sujeto(ismember(etiqueta,clase));                             % seleccion de los trails segun la clase o clases seleccionadas.
        N_trial = length(Izq_Der);                                              % cantidad de trials.
        Canales = [1:size(Izq_Der{1},2)]; N_canal = length(Canales);                            % cantidad de canales.
        temp = spectrogram(Izq_Der{1}(:,1), window, noverlap,Nw);               % tamaño temporal de los datos en cada sujeto.
        temp = size(temp,2);
        X_suj{clase}(su) = {zeros(N_trial,N_canal,floor((Nw/2)+1),temp)};       % Inicializar matricesTrials x canal x frecuencia x tiempo.
        for tri = 1:N_trial
            %             fprintf(['sujeto: ' num2str(su) ' de ' num2str(ns) ' ...trial: ' ...
            %                 num2str(tri) ' de ' num2str(N_trial) '\n'])
            for cnl = 1:N_canal
                signal = Izq_Der{tri}(:,Canales(cnl));                          % signal
                [X_Class, f, t] = spectrogram(signal, window, ...
                    noverlap,Nw,fs);                                            % Calcular STFT
                X_Class = abs(X_Class);                                         % absoluto para tener la señal (freq,tiempo)
                X_suj{clase}{su}(tri,cnl,:,:) = X_Class;                        % Almacenar espectrogramas -> Trials x canal x frecuencia x tiempo
            end
        end
        fprintf(['Sujeto: ' num2str(su) ' de ' num2str(ns) ' ...t: ' ...
            num2str(toc) '\n'])
    end
end
%     % %%%%%%%%%%%%%%%%%%%%%%%%%%% ERD - Normal %%%%%%%%%%%%%%%%%%%%%%%%%%% %% % calculo del ERD
% Reference window
tref = [0,2];
tposref = [find(t==min(t)),find(t==tref(2))];
temp1 = abs(t - t1);                          % ubizacion de tiempo inicial
min1 = min(temp1);
temp2 = abs(t - t2);                          % ubicacion de tiempo final.
min2 = min(temp2);
ul = find(temp1 == min1);                     % ul inicial.
up = find(temp2 == min2);                     % up final.
%
%     r_nc = cell(ns,1);                              % celda para esperanza respecto al tiempo de referencia.
%     r_c = cell(ns,1);                               % celda para esperanza respecto a los trials (intentos).
%     m_c = cell(ns,1);                               % celda para esperanza respecto a los trials (intentos).
%
%     fprintf('ERD in class 1 and 2 \n')
    for i = 1:2                                     % Quantification ERD
        r_nc = cellfun(@(E) squeeze(mean(E(:,:,:,tposref(1):tposref(2)),4)),X_suj{i},...
            'UniformOutput',false);
        r_c = cellfun(@(E) squeeze(mean(E,1)),r_nc,'UniformOutput',false);
        m_c = cellfun(@(E) squeeze(mean(E,1)),X_suj{i},'UniformOutput',false);
        ERD{i} = cellfun(@(A,B) (bsxfun(@times,A,1./B)-1)*100,m_c,r_c,...
            'UniformOutput',false);  % channel x freq x time
    end
%     clear X y
%
%     % Calculate SampEn using a sliding time window of width 2s from 2.5 to 3.5s
%     % for each class.
%     % Analysis window
%     temp1 = abs(t - t3);                            % ubizacion de tiempo inicial
%     min1 = min(temp1);
%     temp2 = abs(t - t4);                            % ubicacion de tiempo final.
%     min2 = min(temp2);
%     ul1 = find(temp1 == min1); ul1 = ul1(1);        % ul inicial.
%     up1 = find(temp2 == min2); up1 = up1(1);        % up final.
%
%     % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SampEn %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
fprintf('SampEn')
w = 2; % seg
w = find(t==w);
% ul1 = find(t == 2.5);
% up1 = find(t == 4.5);
% SampEn1 = Main_SE(ERD{3}{1},ul1,up1,w);       % SampEn Clase 1 sub 1.
% SampEn2 = Main_SE(ERD{3}{2},ul1,up1,w);       % SampEn Clase 2 sub 1.
%
% tt = tes(SampEn1,SampEn2);                    % t-test Clase 1 and 2.
%
% for ch = 1:size(tt,1); ttan(ch) = min(tt(ch,:,2)); end
% P_cor = sort(ttan);
%
SampEn1 = cellfun(@(X) Main_SE(X,ul1,up1,w),ERD{1},'UniformOutput',false);  % ERD class 1
SampEn2 = cellfun(@(X) Main_SE(X,ul1,up1,w),ERD{2},'UniformOutput',false);  % ERD class 2
%
%     fprintf('ttest')
%     tt = cellfun(@(X,Y) tes(X,Y),SampEn1,SampEn2,'UniformOutput',false);        % Ttest - class 1 and 2.
%
%     save('Dataspec.mat','SampEn1','SampEn2','tt','ERD')
%
%     for su = 1:ns
%         for fil = 1:39
%             ttan{su}(fil,:,:) = tt{su}{39};
%         end
%     end
%
%     % load('Data.mat')
%     P_cmin = cellfun(@(X) min(X(:,:,2)),ttan,'UniformOutput',false);             % Pc = min(f)-> Pcf
%     P_cor = cellfun(@(X) sort(X),P_cmin,'UniformOutput',false);                 % Pc_cor = sort(Pc)
%

%     for sub = 1:ns                                                              % relaciona el numero del canal con el orden de los p-valores.
%         for ch = 1:N_canal
%             posch{sub}(ch) = find(P_cmin{sub} == P_cor{sub}(ch));
%         end
%     end
%     save('Datos.mat','P_cor','posch')
% end
%
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%% ERDs-Filter %%%%%%%%%%%%%%%%%%%%%%%%%%%% %% % parametros de la STFT
% if filt == 1
%     ns = numel(X);                                                          % numero de sujetos.
%     N_trials = numel(X{1});
%     [time N_canal] = size(X{1}{1});
%     n_class = 2;
% %     for clase = [1 2]
% %         fprintf(['Clase: ' num2str(clase) '\n'])
% %         for su = 1:ns                                                     % calculo del espectrograma relacionada a EEG.
% %             tic
% %             Sujeto = X{su};                                               % sujeto seleccionado.
% %             etiqueta = y{su};                                             % etiqueta seleccionada.
% %             Dat= Sujeto(ismember(etiqueta,clase));
% %             for k = 1:39
% %                 a = [k k+2];
% %                 Xfreq{k,1} = fcnfiltband(Dat,fs,a,5);                     % funcion Filtro.
% %             end
% %             toc
% %             if clase == 1
% %                 save(['Clase1\Clase1_' num2str(su) '.mat'],'Xfreq','-v7.3')
% %                 clear Xfreq
% %             elseif clase == 2
% %                 save(['Clase2\Clase2_' num2str(su) '.mat'],'Xfreq','-v7.3')
% %                 clear Xfreq
% %             end
% %         end
% %     end
% %
%     % Reference window
% %     tposref = [1,200];
% %     r_nc = cell(ns,1);                              % celda para esperanza respecto al tiempo de referencia.
% %     r_c = cell(ns,1);                               % celda para esperanza respecto a los trials (intentos).
% %     m_c = cell(ns,1);                               % celda para esperanza respecto a los trials (intentos).
%
% %     for clase=[1 2]
% %         fprintf(['ERD in class ' num2str(clase) '\n'])                    % Quantification ERD
% %         for su = 1:ns
% %             tic
% %             load(['Clase' num2str(clase) '\Clase' num2str(clase) '_' num2str(su) '.mat']) % load Xfreq
% %             r_nc = cellfun(@(E) squeeze(mean(E(:,tposref(1):tposref(2),:),2)),Xfreq,'UniformOutput',false);
% %             r_c = cellfun(@(E) squeeze(mean(E,1)),r_nc,'UniformOutput',false);
% %             m_c = cellfun(@(E) squeeze(mean(E,1)),Xfreq,'UniformOutput',false);
% %             ERD{clase}{su} = cellfun(@(A,B) (bsxfun(@times,A,1./B)-1)*100,m_c,r_c,'UniformOutput',false);
% %             clear Xfreq
% %             toc
% %         end
% %     end
% %     clear X y
%    % ----------------------------------------------------------------------
%     load('ERD.mat')
%     % Calculate SampEn using a sliding time window of width 2s from 2.5 to 3.5s
%     % for each class.
%     % Analysis window
%     ul1 = 250;                                        % ul inicial.
%     up1 = 350;                                        % up final.
%
%     % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SampEn %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%     fprintf('SampEn')
%     w = 200; % seg
%     %     w = find(t==w);
%     %     SampEn1 = Main_SE2(ERD{1}{1}{1},ul1,up1,w);                       % SampEn Clase 1 sub 1.
%     %     SampEn2 = Main_SE2(ERD{2}{1}{1},ul1,up1,w);                       % SampEn Clase 2 sub 1.
%     %     tt = tes2(SampEn1,SampEn2);                                       % t-test Clase 1 and 2.
%
%     SampEn1 = cellfun(@(X) cellfun(@(X) Main_SE2(X,ul1,up1,w),X,...
%         'UniformOutput',false),ERD{1},'UniformOutput',false);               % ERD class 1
%     SampEn2 = cellfun(@(X) cellfun(@(X) Main_SE2(X,ul1,up1,w),X,...
%         'UniformOutput',false),ERD{2},'UniformOutput',false);               % ERD class 2
%     fprintf('ttest')
%     tt = cellfun(@(X,Y) cellfun(@(X,Y) tes2(X,Y),X,Y,'UniformOutput',...
%         false),SampEn1,SampEn2,'UniformOutput',false);                      % Ttest - class 1 and 2.
%     save('Datafilter.mat','SampEn1','SampEn2','tt','ERD')
%
%     % -----------------------------------------------------
%     P_cmin = cellfun(@(X) cellfun(@(X) min(X(:,2)),X,'UniformOutput',false),tt,'UniformOutput',false);             % Pc = min(f)-> Pcf
%     P_cor = cellfun(@(X) sort(X),P_cmin,'UniformOutput',false);                 % Pc_cor = sort(Pc)
%
%     for sub = 1:ns                                                              % relaciona el numero del canal con el orden de los p-valores.
%         for ch = 1:N_canal
%             posch{sub}(ch) = find(P_cmin{sub} == P_cor{sub}(ch));
%         end
%     end
%     save('Datos.mat','P_cor','posch')
% end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CSP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
dataset = 0;
load('BCICIV_1');          % raw data
%load('Datos.mat')                                   % data pos.
name = 'BCICIV_1';
fprintf('Loaded data: %s\n',name);
%load('vec.mat')                 % indicadores de los 59 canales
fs = 100;
ns = numel(X);
[~,N_channel] = size(X{1}{1});
sub = 1:ns;
Xd_ = cell(1,ns);
lb = cell(1,ns);
labels = [1 2];
for i = 1:ns
    ind{i} = ismember(y{i},labels);
    Xd_{i} = X{i}(ind{i}); %belonging to labels
    lb{i} = y{i}(ind{i});
end
fprintf(['Clase: ' num2str(labels) ' \n'])
% (MI segment)
tini = 2.5;
tfin = 4.5;
Xdr = fncCutdataf(Xd_,tini,tfin,fs,[8,36]);  % selecciona y filtra

SelecCh = 2;

if SelecCh == 1
    % Analysis
    np = 5; % Kfolds
    [S,Ch,F] = ndgrid(sub,6:59,1:np);
    avtrain = zeros(size(S));
    avtest = zeros(size(S));
    
    for i = 1:numel(S)
        %for sub = 1:numel(X)
        tic
        s = S(i);
        ch = Ch(i);
        f = F(i); % folds
        Xd = Xdr{s};
        tr_ind = cv{s}.training(f); tr_ind = tr_ind(ind{s});
        ts_ind = cv{s}.test(f); ts_ind = ts_ind(ind{s});
        
        Xtrain = Xd(tr_ind);
        Xtest = Xd(ts_ind);
        idtr = find(tr_ind==1);
        ltr = lb{s}(idtr)';
        idts = find(ts_ind==1);
        lts = lb{s}(idts)';
        
        
        
        % CSP
        Xtrain_ = cellfun(@(x) x(:,posch{s}(1:ch)),Xtrain,'UniformOutput',false);
        W_ = fncCSP(Xtrain_,ltr,3);
        
        Xtest_ = cellfun(@(x) x(:,posch{s}(1:ch)),Xtest,'UniformOutput',false);
        
        [Xftrain,Xftest] = fncfeatlog(Xtrain_,Xtest_,W_');
        Mdl = fitcdiscr(Xftrain,ltr); %%%%%%LDA%%%%%%
        predtr = (predict(Mdl,Xftrain));
        predts = (predict(Mdl,Xftest));
        if dataset == 1
            predtr = predtr';
            predts = predts';
        end
        avtrain(i) = sum(predtr == ltr)/numel(ltr);
        avtest(i) = sum(predts == lts)/numel(lts);
        
        %[i,k]
        %end
        %end
        tiempo = toc;
        fprintf('Exp - %d. Time: %f\n',i,tiempo)
        fprintf('Acc    : Training %f. Test %f\n', avtrain(i)*100,avtest(i)*100)
    end
    
    MeanAvtest = mean(avtest,3);
    stdAvtest = squeeze(permute(std(permute(avtest,[3 2 1])),[3 2 1]));
    set(0,'DefaultFigureWindowStyle','docked')
    for suj = 1:ns
        fig=figure(suj);
        errorbar(MeanAvtest(suj,:)',stdAvtest(suj,:)',...
            'CapSize',10,'AlignVertexCenters','on')
        xticks(1:5:59)
        set(fig,'name',['Sujeto' num2str(suj) ': CSP'])
        xticklabels({'6','11','16','21','26','31','36','41','46','51','56'})
        xlabel('N. Channels','Interpreter','latex')
        title(['Subject ' num2str(suj)],'Interpreter','latex')
        xlim([0 55])
        ylim([0.4 1])
        drawnow
    end
    
elseif SelecCh == 2
    
    % Analysis
    np = 5; % Kfolds
    [S,F] = ndgrid(sub,1:np);
    avtrain = zeros(size(S));
    avtest = zeros(size(S));
    
    for i = 1:numel(S)
        %for sub = 1:numel(X)
        tic
        s = S(i);
        %         ch = Ch(i);
        f = F(i); % folds
        Xd = Xdr{s};
        tr_ind = cv{s}.training(f); tr_ind = tr_ind(ind{s});
        ts_ind = cv{s}.test(f); ts_ind = ts_ind(ind{s});
        
        Xtrain = Xd(tr_ind);
        Xtest = Xd(ts_ind);
        idtr = find(tr_ind==1);
        ltr = lb{s}(idtr)';
        idts = find(ts_ind==1);
        lts = lb{s}(idts)';
        
        % CSP
        
        %         Xtrain_ = cellfun(@(x) x(:,vecselec),Xtrain,'UniformOutput',false);
        W = fncCSP(Xtrain,ltr,3);
        
        %         Xtest_ = cellfun(@(x) x(:,vecselec),Xtest,'UniformOutput',false);
        
        [Xftrain,Xftest] = fncfeatlog(Xtrain,Xtest,W');
        Mdl = fitcdiscr(Xftrain,ltr); %%%%%%LDA%%%%%%
        predtr = (predict(Mdl,Xftrain));
        predts = (predict(Mdl,Xftest));
        if dataset == 1
            predtr = predtr';
            predts = predts';
        end
        avtrain(i) = sum(predtr == ltr)/numel(ltr);
        avtest(i) = sum(predts == lts)/numel(lts);
        
        %[i,k]
        %end
        %end
        tiempo = toc;
        fprintf('Exp - %d. Time: %f\n',i,tiempo)
        fprintf('Acc    : Training %f. Test %f\n', avtrain(i)*100,avtest(i)*100)
    end
    MeanAvtest = mean(avtest,2);
    stdAvtest = std(avtest')';
    figure
    errorbar(MeanAvtest',stdAvtest')
    xticks(1:1:7)
    grid on
    xlabel('N. Channels','Interpreter','latex')
    title('CSP','Interpreter','latex')
    drawnow
end

% CSP
% fprintf('Testing Results: %s \n', name)
% for i = 1:numel(sub)
%     acc_s(i,1) = mean(avtest(i,:,:))*100;
%     acc_s(i,2) = std(avtest(i,:,:))*100;
%     fprintf('Acc : %f. +- %f\n', acc_s(i,1),acc_s(i,2))
% end

%         Xs = cellfun(@(x) x(:,posch{sub}(1:ch)),Xd_{sub},'UniformOutput',false);
%
%         %         Xs = cellfun(@(x) squeeze(x(:,band(ii),mask(band(ii),:))),Xb,'UniformOutput',false);
%         C = cellfun(@(x)(cov(x)),Xs,'UniformOutput',false);
%         C = cellfun(@(x)(x/trace(x)),C,'UniformOutput',false);
%         C = cell2mat(reshape(C,[1 1 numel(C)]));
%         W = csp_feats(C(:,:,tr_ind),ys(tr_ind),'train','Q',3);
%
%         tempfeats = csp_feats(C,W,'test');

%A = inv(W);
%feature
%
%
