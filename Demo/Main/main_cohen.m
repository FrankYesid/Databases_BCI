
%% MAIN IWCSP
% COHEN2017: Using spatiotemporal source separation to identify 
% prominent features in multichannel data without sinusoidal filters
%
% L.F. Velasquez-Martinez 
% -------------------------------------------------------------------------

clc, clear %all;
%
% parametros para manejo del codigo
Chan = 8;  % canal (8-C3) - (10-Cz) - (12-C4) se puede colocar cualquiera de los 22.
YYtt = 1;  % 0 - Yt = Wt.Xc, 1 - Yt = Wt*Xc
tam = 500; % 280 = 1seg 500 = 2seg - ventana de analisis.
 
%% Dataset
% dt = 1,2,3:  1:BCICIV_2a
%
% for dt = 1:3
% dt = 1; % load data
% dataset = dt;

% dataset =0;
% load('G:\work_dic\Articulo2016\codes\BCICIIIIVa.mat');
% name = 'BCICIIIVa';
% fs = 100;
% sub = 1:5; %% Subject to analyze
% Xd_ = cell(1,length(sub));
% lb = cell(1,length(sub));
% labels = [1 2];
% for i = 1:length(sub)
%     ind{i} = ismember(y{i},labels);
%     Xd_{i} = X{i}(ind{i});
%     lb{i} = y{i}(ind{i});
% end
% % (MI segment)
% tini = 0.0001;
% tfin = 3.5;
% Xdr = fncCutdataf(Xd_,tini,tfin,fs,[8,30]);  % ventana de tiempo
% fprintf('Loaded data: %s\n',name);
% ns = size(Xd_,2);

%%
load('H:\ERD\BCICIV_2a.mat');% raw data
name = 'BCICIV2a';
dataset = 1;
ns = size(X,1);
sub = 1:ns;
Xd_ = cell(1,ns);
lb = cell(1,ns);
labels = [1 2]; % clases
for i = 1:ns
    ind{i} = ismember(y{i},labels);
    Xd_{i} = X{i}(ind{i}); %belonging to labels
    lb{i} = y{i}(ind{i});
end
fprintf('Loaded data: %s\n',name);

% %% database - EEGLAB
% da = cell(1,1);
% s = 5;
% for j = 1:numel(Xd_{s})
%     da{1}(:,:,j) = Xd_{s}{j}';
% end
% da=permute(da,[3 2 1]);
% save(['BCICIV_2a' num2str(s) 'class1&2.mat'],'da')

%% Filter Band
clc;
a = [8 30];                                  % rango del filtro.
A = cell(9,1);                               % celdas de la señal filtrada.
filt = 0;
% filtrado de cad sujeto.
if filt == 1
    for k = 1:9
        fprintf(['Sujeto: ' num2str(k) ' de 9' '\n'])
        Xfreq = fcnfiltband(X{k},fs,a,5);           % funcion Filtro.
        A{k} = Xfreq;
    end
    X = A;
end

%% Channels C3, C4                                                         
Xa_ = cell(1,ns);
for i = 1:ns
   Xa_{i} = cellfun(@(X) X(:,Chan),Xd_{i},'UniformOutput',false);
end
Xd_ = Xa_; clear Xa_;


%% Selecting data for analysis
% (MI segment)
tini = 2;                             % tiempo de señal
if tam == 250
    trfin = 3;
else 
    trfin = 4;
end
samples = (trfin-tini)*fs;
Xdd = fncCutdata(Xd_,tini,trfin,fs);   % selecciona y filtra 1x500
for s = 1: ns
    Xdc{s} = cellfun(@(X) X-mean(X),Xdd{s},'UniformOutput',false);
end
% (Reference segment)
trini = 0;                            % tiempo de referencia.
if tam == 250
    trfin = 1;
else 
    trfin = 2;
end
samplesr = (trfin-trini)*fs;
Xdr = fncCutdata(Xd_,trini,trfin,fs);  % selecciona y filtra 1x500
for s = 1: ns
    Xrc{s} = cellfun(@(X) X-mean(X),Xdr{s},'UniformOutput',false);
end

%% hankelizar - time. 1x500
datd = fncHank(Xdc); % data
datr = fncHank(Xrc); %reference
Xdd = datd; clear datd;
Xdr = datr; clear datr;

%% mean


%% Regularization
for i = 1:ns
   datd{i} = cellfun(@(X) X + 0.1.*eye(tam).*X,Xdd{i},'UniformOutput',false);
   datr{i} = cellfun(@(X) X + 0.1.*eye(tam).*X,Xdr{i},'UniformOutput',false);
end
Xdd = datd; clear datd;
Xdr = datr; clear datr;

%% hankelizar todo el tiempo. 1x1750 
% datd = cell(ns,1);
% for y = 1:ns
%    datd{y} = cell(numel(Xd_{y}),1);
%    for i = 1:numel(Xd_{y})
%        datd{y,1}{i} = Xd_{y}{i};
%    end
% end
% Xd_ = datd; clear datd;
% DatXc3 = fncHank2(Xd_); % data - time serial

%% Analysis
np = 5; % Kfolds
[S,F] = ndgrid(sub,1:np);
avtrain = zeros(size(S));
avtest = zeros(size(S));

% pruebas = numel(S);

for i = 1: numel(S)
%     tic
    s = S(i);
    f = F(i); % fold
    
    Xd = Xdd{s};%Xdata  ----
    Xr = Xdr{s};%Xreference
    
    tr_ind = cv{s}.training(f); tr_ind = tr_ind(ind{s});
    ts_ind = cv{s}.test(f); ts_ind = ts_ind(ind{s});
    
    Xtrd = Xd(tr_ind);
    Xtrr = Xr(tr_ind);
    
    Xtest = Xd(ts_ind);
    idtr = find(tr_ind==1);
    ltr = lb{s}(idtr)';
    idts = find(ts_ind==1);
    lts = lb{s}(idts)';
    
    
    %% Regularization
%     datd = cellfun(@(X) X+(0.1*eye(tam).*X),Xtrd,'UniformOutput',false);
%     datr = cellfun(@(X) X+(0.1*eye(tam).*X),Xtrr,'UniformOutput',false);
%     Xtrd = datd; clear datd;
%     Xtrr = datr; clear datr;
    
    %%
    %Etapa de diseño
    Wt = fncCSPcohen(Xtrd,Xtrr); %Finding temporal filters 
    %[Xd_,Xr_] = fncrotate(Xtrd,Xtrr,Ws); %Rotating data
    %[Xhd_,Xhr_] = fncHankel(Xd_,Xr_,500);
    %Wt = fncCSPcohen(Xhd_,Xhr_); %Finding temporal filters 
    
    % plot Wt
    figure,plot(Wt),title('Wt')
    % plot fft de Wt.
    figure,plot(0:1023,smooth(abs(fft(Wt,1024)),15,'loess')),xlim([0 60]),title('FFT ')
    
    %% selecciona multiplicacion o convolucion en la Yt
    if YYtt == 0
        %% Yt = Wt.Xc
        Yt = cell(ns,1);
        fprintf('Yt = Wt.Xc \n')
        for i = 1:ns
            load(['dd' num2str(i)  '.mat'])
            Yt{i} = cell(1,numel(dd{i}));
            fprintf(['Sujeto: ' num2str(i) ' de 9 \n'])
            for k = 1:numel(dd{i})
                Yt{i}{k} = cell(1,7);
                ja=1;
                for j = 1:250:1501
                    Yt{i}{k}{ja} = Wt*dd{i}{k}(j:j+249,1:250); % multiplica cada una de las ventanas de 1seg o 2 seg 
                    ja=ja+1;
                end
            end
            clear dd
        end
        % guardar matriz 
        save('Yt','Yt')
        
        %% Mean - std
        load('Yt.mat')
        ytm = cell(ns,1);
        fprintf('Mean - Std \n')
        for i = 1:ns
            fprintf(['Sujeto: ' num2str(i) ' de 9 \n'])
            for j = 1:numel(Yt{i})
                ytm{i}(j,:,:) = cell2mat(Yt{i}{j}'); % agrupa los diferentes intentos (trials) en una sola matriz
            end
            for l = 1:7
                channel = ytm{i}(:,l,:);
                mYt{i,1}(l,:) = mean(channel); % calcula la media por todos los trials 
                sYt{i,1}(l,:) = std(channel);  % calcula la 
            end
        end
        % guardar matrices 
        save('ytm.mat','ytm')
        save('mYt.mat','mYt')
        save('sYt.mat','sYt')
        
        %% Image Mean - STD todos los sujetos
        set(0,'DefaultFigureWindowStyle','docked')
        for s=1:ns                                  
            if Chan == 8
                nam = 'C3';
            else
                nam = 'C4';
            end
            load('mYt.mat')
            load('sYt.mat')
            da = 1;
            t=0:1/250:7-1/250;
            for i = 1:250:1501
                MYT(i:i+249) = mYt{s}(da,:);
                SYT(i:i+249) = sYt{s}(da,:);
                da = da+1;
            end
            figure1 = figure('name','Errobar');
            hold on
            errorbar(t,MYT,SYT)
            title(['Sujeto: ' num2str(s) ' Canal: ' nam])
            
            plot(t.*0+2,t.*100,'r--','LineWidth',2.5)
            plot(t.*0+3.25,t.*100,'b--','LineWidth',2.5)
            plot(t.*0+6,t.*100,'m--','LineWidth',2.5)
            
            plot(t.*0+2,-t.*100,'r--','LineWidth',2.5)
            plot(t.*0+3.25,-t.*100,'b--','LineWidth',2.5)
            plot(t.*0+6,-t.*100,'m--','LineWidth',2.5)
            ylim([-500 500])
            hold off
            
            figure2 = figure('name','Mean');
            hold on
            plot(t',MYT)
            plot(t.*0+2,t.*100,'r--','LineWidth',2.5)
            plot(t.*0+3.25,t.*100,'b--','LineWidth',2.5)
            plot(t.*0+6,t.*100,'m--','LineWidth',2.5)
            
            plot(t.*0+2,-t.*100,'r--','LineWidth',2.5)
            plot(t.*0+3.25,-t.*100,'b--','LineWidth',2.5)
            plot(t.*0+6,-t.*100,'m--','LineWidth',2.5)
            ylim([-500 500])
            hold off
            title(['Sujeto: ' num2str(s) ' Canal: ' nam])
        end
    else
        %% Yt = Wt*Xc
        Yt = cell(ns,1);
        fprintf('Yt = Wt*Xc')
        for i = 1:ns
            fprintf(['Sujeto: ' num2str(i) ' de 9 \n'])
            Yt{i} = cellfun(@(X) conv(X,Wt','same')',Xd_{i},'UniformOutput',false);
        end
        % guardar matriz
        save('Yt','Yt')
        
        %% Mean - STD 
        load('Yt.mat')
        fprintf('Mean - Std \n')
        for i = 1:ns
            fprintf(['Sujeto: ' num2str(i) ' de 9 \n'])
            for j = 1:numel(Yt{i})
                mYt{i} = mean(cell2mat(Yt{i}),1);
                sYt{i} = std(cell2mat(Yt{i}),1);
            end
        end
        % guardar matriz de mean - std
        save('mYt.mat','mYt')
        save('sYt.mat','sYt')
        
        %% Image  Mean - STD todos los sujetos
        set(0,'DefaultFigureWindowStyle','docked')
        for s=1:ns                                   % seleccionar sujeto
            if Chan == 8
                nam = 'C3';
            else
                nam = 'C4';
            end
            load('mYt.mat')
            load('sYt.mat')
            
            t=0:1/250:7-1/250;
            MYT = mYt{s};
            SYT = sYt{s};
                
            figure1 = figure('name','Errobar');
            hold on
            errorbar(t,MYT,SYT)
            title(['Sujeto: ' num2str(s) ' Canal: ' nam])
            
            plot(t.*0+2,t.*100,'r--','LineWidth',2.5)
            plot(t.*0+3.25,t.*100,'b--','LineWidth',2.5)
            plot(t.*0+6,t.*100,'m--','LineWidth',2.5)
            
            plot(t.*0+2,-t.*100,'r--','LineWidth',2.5)
            plot(t.*0+3.25,-t.*100,'b--','LineWidth',2.5)
            plot(t.*0+6,-t.*100,'m--','LineWidth',2.5)
            ylim([-500 500])
            hold off
            
            figure2 = figure('name','Mean');
            hold on
            plot(t',MYT)
            plot(t.*0+2,t.*100,'r--','LineWidth',2.5)
            plot(t.*0+3.25,t.*100,'b--','LineWidth',2.5)
            plot(t.*0+6,t.*100,'m--','LineWidth',2.5)
            
            plot(t.*0+2,-t.*100,'r--','LineWidth',2.5)
            plot(t.*0+3.25,-t.*100,'b--','LineWidth',2.5)
            plot(t.*0+6,-t.*100,'m--','LineWidth',2.5)
            ylim([-250 250])
            hold off
            title(['Sujeto: ' num2str(s) ' Canal: ' nam])
        end
    end
    
%     %%
%     %Aplicando rotaciones sobre Xtrain y Xtest
%     [Xdtrain,Xdtest] = fncrotate(Xtrd,Xtest,Wt); %Rotating data
%     [Xdftrain,Xdftest] = fncKernelfiltercohen(Xdtrain,Xdtest,Wt); %feature
%     
%     Mdl = fitcdiscr(Xdftrain,ltr); %%%%%%LDA%%%%%%
%     predtr = (predict(Mdl,Xdftrain));
%     predts = (predict(Mdl,Xdftest));
%     
%     if dataset == 1
%         predtr = predtr';
%         predts = predts';
%     end
%     avtrain(i) = sum(predtr == ltr)/numel(ltr);
%     avtest(i) = sum(predts == lts)/numel(lts);
%     tiempo = toc;
%     fprintf('Exp - %d. Time: %f\n',i,tiempo)
%     fprintf('Acc    : Training %f. Test %f\n', avtrain(i)*100,avtest(i)*100)
%     %[i,k]
end
% % % CSP
% fprintf('Testing Results: %s \n', name)
% for i = 1:numel(sub)
%     acc_s(i,1) = mean(avtest(i,:,:))*100;
%     acc_s(i,2) = std(avtest(i,:,:))*100;
%     fprintf('Acc : %f. +- %f\n', acc_s(i,1),acc_s(i,2))
% end

% 

%% 
% load 

%%
% %%%%%%%%%%%%%%%%%%%%%% Graphs

% save('results_CSP1.mat','av_s','std_s');
% figure;
% errorbar(1:9,av_s,std_s,'-s','LineStyle','none','MarkerSize',10,...
%     'MarkerEdgeColor','blue','MarkerFaceColor','blue'), ylim([0 1])
% figure
% bar(1:9,av_s)
% ylim([0 1])
%
% %% TRCSP (alpha)
% for i = 1:numel(sub)
%     for j = 1:numel(alph)
%     av_s(j,i) = mean(mean(squeeze(avtest(i,:,j))));
%     std_s(j,i) = std((squeeze(avtest(i,:,j))));
%     end
% end
% save('results_WTRCSP.mat','av_s','std_s','weightVector');
% figure;
% for j = 1:9
%     subplot(3,3,j)
%     fncGraph_acc(av_s(:,j)',std_s(:,j)',alph,0,1);
%     title(['S',num2str(j)]);
% end
% %% WTRCSP
% figure;
% for j = 1:9
%     plot(1:9,weightVector{j},'-s','MarkerSize',10,...
%      'MarkerEdgeColor','blue','MarkerFaceColor','blue'); ylim([0 1]);
%     title(['S',num2str(j)]);
% end
