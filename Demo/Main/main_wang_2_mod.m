clear all; close all; clc
%%

% % Direccion de la base de datos
% SUBJECTS_DIR = 'D:\BCI';
SUBJECTS_DIR = 'F:\BCI';
SUBJECTS_DIR2 = 'G:\Dropbox\ERD\results_ERDfc_subjects\BCI';
% % Direccion del fold de las funciones
addpath(genpath('C:\Users\lfvelasquezm\Dropbox\ERD\Codes\TP\Matlab_wang\csp\CSP_fun\functions'))
% addpath(genpath('C:\Users\lfvelasquezm\Desktop\frank\functions'))
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
SS = [15 37 32 12 18 42 34 3 7 35 33 21 2 4 39 29 43 28]; % UNO BUENO Y UNO MALO%%%%INDEXACDOS DE ACIERDO A CSP
% SS = [37 32 12 18 42 34 3 7];
% if strcmp(COHORT,'GIGASCIENCE_')
%     SubInd = [50,14];
%     SS(SubInd) = [];
% end

%% paramaters definition
% Rep = 10000;  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% numero de repeticiones
% PPval = cell(numel(SS),1);
% rho = [];
% tstart = 0;
% tend = 2;
w = 500;            % tamaño de la ventana.

% definir parametros de filter bank
f_low  = 4;
f_high = 40;
Window = 2;
Ovrlap = 1;
filter_bank = [f_low:Ovrlap:f_high-Window;...
    f_low+Window:Ovrlap:f_high]';

%%
for s = SS
    %     clearvars -except s SS rho experiment_name COHORT param SUBJECTS SUBJECTS_DIR Acc table PPval Rep tstart tend
    load([SUBJECTS_DIR filesep SUBJECTS{s} filesep 'eeg' filesep 'raw.mat'])
    labels = unique(y);
    N_channel = size(X{1},2);
    N_class = numel(labels);
    ind = ismember(y,labels);
    y = y(ind);
    X = X(ind);
    X = cellfun(@(x) double(x) ,X,'UniformOutput',false);
    tic
    for clase = 1:N_class
        fprintf(['Clase: ' num2str(clase) '\n'])
        %         Dat = X(ismember(y,clase));
        temp = cellfun(@(x) x.^2, X(ismember(y,clase)),'UniformOutput',false);
        for b = 1:size(filter_bank,1)
            Xfreq{b,clase} = fcnfiltband_vector(temp,fs,filter_bank(b,:),5);   % funcion Filtro.
        end
    end
    
    N_filt = size(Xfreq,1); % Numero de filtros.
    r_nc = cell(N_filt,1);  % celda para esperanza respecto al tiempo de referencia.
    r_c = cell(N_filt,1);   % celda para esperanza respecto a los trials (intentos).
    m_c = cell(N_filt,1);   % celda para esperanza respecto a los trials (intentos).
    ERD = cell(1,N_class);  % ERDs por clase.
    tposref = [1,fs*2];     % tiempo de referencia para el calculo del ERDs.
    ul1 = seg_start;        % inicial de la posicion de la ventana deslizante.
    up1 = seg_end;          % final de la posicion de la ventana deslizante.
    
    for clase = 1:N_class
        fprintf(['ERD in class ' num2str(clase) '\n']) % Quantification ERD
        r_nc = cellfun(@(E) squeeze(mean(E(:,tposref(1):tposref(2),:),2)),Xfreq(:,clase),'UniformOutput',false);
        r_c = cellfun(@(E) squeeze(mean(E,1)),r_nc,'UniformOutput',false);
        m_c = cellfun(@(E) squeeze(mean(E,1)),Xfreq(:,clase),'UniformOutput',false);
        ERD{clase} = cellfun(@(A,B) (bsxfun(@times,A,1./B)-1)*100,m_c,r_c,'UniformOutput',false);
    end
    fprintf(['ERD done Sub: ' num2str(s) '\n'])
    fprintf(['Sample Entropy Sub: ' num2str(s) '\n'])
    SampEnc1{s} = Main_SE2(ERD{1},ul1,up1,w);  % ERD class 1
    SampEnc2{s} = Main_SE2(ERD{2},ul1,up1,w);  % ERD class 2
    fprintf(['Sample Entropy done Sub: ' num2str(s) '\n'])
    
    %     fprintf('Sub - %d. Fold - %d. Time: %f\n',s,f,time/60)
    
    fprintf(['T-test Sub: ' num2str(s) '\n'])
    tt{s} = tes2(SampEnc1{s},SampEnc2{s});
    fprintf(['T-test done Sub: ' num2str(s) '\n'])
    
    P_cmin{s} =  min(tt{s});                % Pc = min(f)-> Pcf
    P_cor{s} = sort(P_cmin{s});             % Ordenar de menor a mayor.
    
    for ch = 1:N_channel
        temp = find(P_cmin{s} == P_cor{s}(ch));
        if length(temp) > 1
            temp = temp(1);
            posch{s}(ch) = temp;
            P_cmin{s}(temp)= 1000;
        else
            posch{s}(ch) = temp;
        end
    end
    
    % Guardar resultados
    save([SUBJECTS_DIR filesep SUBJECTS{s} filesep 'results\' experiment_name 'posch.mat'],'posch','SampEnc1','SampEnc2','tt','P_cmin');
    % %% Cargar mascara
    %     load([SUBJECTS_DIR filesep SUBJECTS{s} filesep 'results\' experiment_name 'rho.mat'],'rho_s');
    %     if sum(isnan(rho_s))>1; rho_s(isnan(rho_s)==1)=1;  end %para los nan detectados
    %     threshold = linspace(min(rho_s),max(rho_s),100);
    
    %     X = cellfun(@(x) double(x) ,X(ind),'UniformOutput',false);
    %     clc
    ind = ismember(y,labels);
    acc=nan(5,N_channel,numel(param));
    ks=nan(5,numel(threshold),numel(param));
    Xcp = cell(5,1);
    for fold = 1:5
        %         tic;
        tr_ind   = cv.training(fold); tr_ind = tr_ind(ind);
        ts_ind   = cv.test(fold); ts_ind = ts_ind(ind);
        %         for u = 1: numel(threshold)
        tic
        %             mask = reshape(p,[size(Pval,1) size(Pval,2)]) < threshold(u);
        %             mask = reshape(p,[33,22])<= threshold(u);
        %             SelBand = sum(mask,2);
        %             valQ = floor(SelBand/2);
        %             selected = SelBand>=2;
        %             band = find(selected==1);
        %             if sum(selected) == 0
        %                 continue
        %             end
        for ch = 6:N_channel
            %             [~,chan]=find(mask(band(b),:)>0);
            %             if numel(chan) <= 1
            %                 continue
            %             end
            Xc = [];
            Xa = cellfun(@(x) x(:,posch(1:ch)),X,'UniformOutput',false);
            %             temp = fcnfiltband(Xa, fs, filter_bank(b,:), 5);
            Xb = cellfun(@(x) x(seg_start:seg_end,:),Xa','UniformOutput',false);
            C =  cell2mat(reshape(cellfun(@(x)(cov(x)/trace(cov(x))),Xb,'UniformOutput',false),[1 1 numel(Xb)]));
            W = csp_feats(C(:,:,tr_ind),y(tr_ind),'train','Q',floor(numel(chan)/2));
            Xc = [Xc  csp_feats(C,W,'test')];            
            clear Xb C
            Xcp{fold} = Xc;
            % Lasso
%             target = mapminmax(y(tr_ind)')';
%             B = lasso(Xc(tr_ind,:),target,'Lambda',param);
%             selected_feats = abs(B)>eps;
%             for l=1:numel(param)
%                 Xcc{l} =  Xc(:,selected_feats(:,l));
%             end
%             clear Xc
            
%             for l=1:numel(param)
%                 if size(Xcc{l},2)<2
%                     continue
%                 end
                mdl = fitcdiscr(Xc(tr_ind,:),y(tr_ind));
                acc(fold,ch) = mean(mdl.predict(Xc(ts_ind,:))==reshape(y(ts_ind),[sum(ts_ind) 1]));
                %Confusion Matrix
                tar_pred = mdl.predict(Xc(ts_ind,:)); %tar_pred(tar_pred==1)=0; tar_pred(tar_pred==2)=1;
                tar_true = reshape(y(ts_ind),[sum(ts_ind) 1]); %tar_true(tar_true==1)=0; tar_true(tar_true==2)=1;
                conM = confusionmat(tar_true,tar_pred);
                ks(fold,ch) = kappa(conM);
                %plotconfusion(tar_true',tar_pred');
%             end %lambda
        end % channels
        clear Xc
        fprintf(['Threshold...' num2str(u) '...' num2str(toc) '\n'])
        %             [fold,u]
        %         end % signficance
        fprintf(['Sujeto: ' num2str(s) ' Fold...' num2str(fold) '\n'])
        toc
    end % folds
     act = squeeze(mean(acc,1));
    [dato,nch_s] = max(act);
    actstd = std(acc,1); actstd = actstd(nch_s);
    table(s,:) = [nch_s,dato*100,actstd*100];
    toc
    %% Guardar resultados
    save([SUBJECTS_DIR filesep SUBJECTS{s} filesep 'results' filesep experiment_name 'Results.mat'],'acc','table','ks','Xcp');
    save([SUBJECTS_DIR2 filesep SUBJECTS{s} filesep 'results' filesep experiment_name 'Results.mat'],'acc','table','ks');
    fprintf([' ...acc: ' num2str(mmacc(s)*100,'%02.1f') ...
        ' ...time: ' num2str(toc) '\n']);
end
diary('')