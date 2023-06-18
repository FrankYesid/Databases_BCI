clear ; clc
% load('BCICIV_1.mat')   % carga base de datos EEG. BCICIV_1
load('BCICIV_2a.mat')
name  = 'BCIIV_2a';
fprintf('Loaded data: %s\n',name);
% load('vec.mat')                 % indicadores de los 59 canales
% fs = 100;
ns = numel(X);
[~,N_channel] = size(X{1}{1});
sub = 1:9;
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
% Xdr = fncCutdataf(Xd_,tini,tfin,fs,[2,40]);  % selecciona y filtra
Xdr = Xd_;

% filter - bank
bw = 2;
freq = (4+bw:2:40-bw)';
filter_bank = [freq-bw freq+bw];

% Analysis
np = 5; % Kfolds
[S,F] = ndgrid(sub,1:np);
avtrain = cell(size(S));
avtest = cell(size(S));

SampEnc1 = cell(ns,np);
SampEnc2 = cell(ns,np);
% tt = cell(ns,1);
% P_cmin = cell(ns,1);
% P_cor = cell(ns,1);
% posch = cell(ns,np);
dkl_w = cell(ns,np);
jray_w = cell(ns,np);
% clc
% load('tt_final.mat')
% P_cmintemp = cellfun(@(X) min(X),tt,'UniformOutput',false);  % Pc = min(f)-> Pcf
% P_cmin = cell(ns,1);
% for s = 1:ns
%     P_cmin{s} = cell2mat(P_cmintemp(s,:)');
% end
%
% P_cmin = cellfun(@(X) mean(X),P_cmin,'UniformOutput',false);
%
% P_cor = cellfun(@(x) sort(x),P_cmin,'UniformOutput',false);   %ordena de menor a mayor
%
% % relaciona el numero del canal con el orden de los p-valores.
% for s = 1:ns
%     for ch = 1:N_channel
%         tmp = find(P_cmin{s} == P_cor{s}(ch));
%         if length(tmp)>1
%             tmp = tmp(1);
%             P_cmin{s}(tmp) = 1000;
%         end
%         posch{s}(ch) = tmp;
%     end
% end
% load('ResultsSamEnBCIIV1_wang.mat')
for i = 1:numel(S)
    s = S(i);
    f = F(i); % folds
    Xd = Xdr{s};
    tic
    N_class = numel(labels);
    tr_ind = cv{s}.training(f); tr_ind = tr_ind(ind{s});
    ts_ind = cv{s}.test(f); ts_ind = ts_ind(ind{s});
    tic
    Xtrain = Xd(tr_ind);
    N_train = numel(Xtrain);
    Xtest = Xd(ts_ind);
    idtr = find(tr_ind==1);
    ltr = lb{s}(idtr)';
    idts = find(ts_ind==1);
    lts = lb{s}(idts)';
    for clase = 1:N_class
        fprintf(['Clase: ' num2str(clase) '\n'])
        Sujeto = Xtrain;  %sujeto seleccionado.
        etiqueta = ltr;   % etiqueta seleccionada.
        Dat = Sujeto(ismember(etiqueta,clase));
        temp = cellfun(@(x) x.^2, Dat,'UniformOutput',false);
        parfor b = 1:numel(freq)
            Xfreq{b,clase} = fcnfiltband_vector(temp,fs,filter_bank(b,:),5);   % funcion Filtro.
        end
    end
    %
    N_filt = size(Xfreq,1);
    tposref = [1,200];                                   % Reference window
    r_nc = cell(N_filt,1); % celda para esperanza respecto al tiempo de referencia.
    r_c = cell(N_filt,1);  % celda para esperanza respecto a los trials (intentos).
    m_c = cell(N_filt,1);  % celda para esperanza respecto a los trials (intentos).
    ERD = cell(1,N_class);
    toc
    tic
    parfor clase=1:2
        fprintf(['ERD in class ' num2str(clase) '\n']) % Quantification ERD
        
        r_nc = cellfun(@(E) squeeze(mean(E(:,tposref(1):tposref(2),:),2)),Xfreq(:,clase),'UniformOutput',false);
        r_c = cellfun(@(E) squeeze(mean(E,1)),r_nc,'UniformOutput',false);
        m_c = cellfun(@(E) squeeze(mean(E,1)),Xfreq(:,clase),'UniformOutput',false);
        ERD{clase} = cellfun(@(A,B) (bsxfun(@times,A,1./B)-1)*100,m_c,r_c,'UniformOutput',false);
    end
    toc
    % Calculate SampEn using a sliding time window of width 2s from 2.5 to 3.5s for each class.
    ul1 = 250;  % inicial de la posicion de la ventana deslizante.
    up1 = 350;  % final de la posicion de la ventana deslizante.
    w = 200;    % seg
    
    fprintf('SampEn \n')
    SampEnc1{s,f} = Main_SE2(ERD{1},ul1,up1,w);  % ERD class 1
    SampEnc2{s,f} = Main_SE2(ERD{2},ul1,up1,w);  % ERD class 2
    %     time = toc;
    fprintf(['Sample Entropy done Sub: ' num2str(s) ' Fold:' num2str(f) '\n'])
    %     fprintf('Sub - %d. Fold - %d. Time: %f\n',s,f,time/60)
    
    fprintf('ttest \n')
    tt{s,f} = tes2(SampEnc1{s,f},SampEnc2{s,f});
    
    P_cmin{s,f} =  min(tt{s,f});    % Pc = min(f)-> Pcf
    P_cor{s,f} = sort(P_cmin{s,f});     % Ordenar de menor a mayor.
    
    for ch = 1:N_channel
        temp = find(P_cmin{s,f} == P_cor{s,f}(ch));
        if length(temp) > 1
            temp = temp(1);
            posch{s,f}(ch) = temp;
            P_cmin{s,f}(temp)= 1000;
        else
            posch{s,f}(ch) = temp;
        end
    end
    save('ResultsSamEnBCIIV2a_wang.mat','posch','SampEnc1','SampEnc2','tt','P_cmin')
    
    for jj = 6:N_channel % CSP
        tic
        Xtrain_ = cellfun(@(x) x(:,posch{s,f}(1:jj)),Xtrain,'UniformOutput',false);
        W_ = fncCSP(Xtrain_,ltr,3);
        
        Xtest_ = cellfun(@(x) x(:,posch{s,f}(1:jj)),Xtest,'UniformOutput',false);
        
        [Xftrain,Xftest] = fncfeatlog(Xtrain_',Xtest_',W_');
        Mdl = fitcdiscr(Xftrain,ltr); %%%%%%LDA%%%%%%
        predtr = (predict(Mdl,Xftrain));
        predts = (predict(Mdl,Xftest));
        
        avtrain{i}(jj) = sum(predtr == ltr')/numel(ltr);
        avtest{i}(jj) = sum(predts == lts')/numel(lts);
        
        tiempo = toc;
        fprintf('Suj: %d. Fold: %d. Relevence ch - %d.   Time: %f\n',s,f,jj,tiempo/60)
        fprintf('Acc    : Training %f. Test %f\n', avtrain{i}(jj)*100,avtest{i}(jj)*100)
    end
    save('ResultsAccBCICIV_wang.mat','avtest')
    
    %     for ii = 6:59
    %         tic
    %         Xtrain_ = cellfun(@(x) x(:,posch{s}(1:ii)),Xtrain,'UniformOutput',false);
    %         Xtest_ = cellfun(@(x) x(:,posch{s}(1:ii)),Xtest,'UniformOutput',false);
    %         tini = 2.5; % (MI segment)
    %         trfin = 4.5;
    %         samples = (trfin-tini)*fs;
    %         Xtrain_ = fncCutdataf(Xtrain_,tini,trfin,fs,[8,30]);
    %         Xtest_ = fncCutdataf(Xtest_,tini,trfin,fs,[8,30]);
    %         %         Xtrain_ = fncCutdata_(Xtrain_,tini,trfin,fs);   % selecciona
    %         %         Xtest_ = fncCutdata_(Xtest_,tini,trfin,fs);   % selecciona
    %
    %         %         Xtrain_ = cellfun(@(Xa) Xa-mean(Xa),Xtrain_,'UniformOutput',false); % media cero
    %         %         Xtest_ = cellfun(@(Xa) Xa-mean(Xa),Xtest_,'UniformOutput',false); % media cero
    %
    %         [W,Cov] = fncCSP(Xtrain_,ltr,3);
    %
    %         %         parfor ii = 1:numel(F) %frequencies and time
    %         %         Xd_ = X{ii};
    %         %         [W,Cov] = fncCSP(Xd_,ytr,3); %------------ W rotation CSP-based
    %         %------------ Rayleight quotient using the first component
    %         %         jtemp(ii) = trace(W*Cov{1}*W')/trace(W*Cov{2}*W');
    %
    %         %feature
    %         [Xftrain,Xftest,covXtest] = fncfeatlog(Xtrain_,Xtest_,W');
    %         Xf1 = Xftest(lts==1,:);
    %         Xf2 = Xftest(lts==2,:);
    %         Cov{1} = covXtest(:,:,lts==1); Cov{1}=mean(Cov{1},3);
    %         Cov{2} = covXtest(:,:,lts==2);  Cov{2}=mean(Cov{2},3);
    %
    %         %Div
    %         k = 6;
    %         %for ff = 1 :6
    %         tmp1 = trace((cov(Xf1))\cov((Xf2)));
    %         tmp2 = (mean(Xf1)-mean(Xf2))*inv(cov(Xf2))*(mean(Xf1)-mean(Xf2))';
    %         tmp3 = log(det(cov(Xf1))/det(cov(Xf2)));
    %         tdkl = 0.5*(tmp1+tmp2-k+tmp3);
    %         %end
    %         %dkl{jj}(ii) = sum(tdkl);
    %         dkl_w{i}(ii) = tdkl;
    %
    %         %Rayleight
    %         %         j(ii) = (W(1,:)*Cov{1}*W(1,:)')/(W(1,:)*Cov{2}*W(1,:)');
    %         jtemp(ii) = (W(:,1)'*Cov{1}*W(:,1))/(W(:,1)'*Cov{2}*W(:,1));
    %         %trace(W*Cov{1}*W')/trace(W*Cov{2}*W');
    %         Mdl = fitcdiscr(Xftrain,ltr); %%%%%%LDA%%%%%%
    %         predtr = (predict(Mdl,Xftrain));
    %         predts = (predict(Mdl,Xftest));
    %
    %         avtrain{i}(ii) = sum(predtr == ltr)/numel(ltr);
    %         avtest{i}(ii) = sum(predts == lts)/numel(lts);
    %         time = toc;
    %         fprintf('Suj: %d. Fold: %d. Relevence ch - %d.   Time: %f\n',s,f,ii,time)
    %         fprintf('Acc    : Training %f. Test %f\n', avtrain{i}(ii)*100,avtest{i}(ii)*100)
    %     end
    %     jray_w{i} = jtemp;
    %     save('ResultsAccBCIIV1_f2.mat','avtest','dkl_w','jray_w')
    
end

% load('ResultsAccBCIIV1_f2.mat')
% set(0,'DefaultFigureWindowStyle','docked')
% a=1;
% for s = [1,2,6,7] %[3,4]
%     for i =1:5; Dkl{s}(i,:) = dkl_w{s,i}/max(dkl_w{s,i}); end
%     for i =1:5; j_{s}(i,:) = jray_w{s,i}/max(jray_w{s,i}); end
%     for i =1:5; ac{s}(i,:) = avtest{s,i}; end
%     figure(a)
%     a = a+1;
%     subplot(3,1,1)
%     errorbar(mean(ac{s}),std(ac{s})); title(['Acc' ' subject: ' num2str(s)])
%     ylim([0 1])
%     subplot(3,1,2)
%     errorbar(mean(Dkl{s}),std(Dkl{s})); title(['Div' ' subject: ' num2str(s)])
%     ylim([0 1])
%     subplot(3,1,3)
%     errorbar(mean(j_{s}),std(j_{s})); title(['Rayleight' ' subject: ' num2str(s)])
%     ylim([0 1])
% end

% hold on
% subplot(3,1,1)
% errorbar(mean(ac{s}),std(ac{s}),'r'); title(['Acc' ' subject: ' num2str(s)])
%
% cc=1;
% for s = [1,2,6,7]; %[max(mean(Dkl_w{s})),std(Dkl{s}(:,find(mean(Dkl_w{s})==max(mean(Dkl_w{s}))))),find(mean(Dkl_w{s})==max(mean(Dkl_w{s})))]
%     [M(cc,1),Ich(cc,1)] = max(mean(ac{s}));
%     M(cc,2) = std(ac{s}(:,Ich(cc,1)));
%     [max(mean(ac{s})),std(ac{s}(:,Ich(cc,1))),find(mean(ac{s})==max(mean(ac{s})))];
%     cc= cc+1;
% end
%
%
% % type 1
% %numero de canales
%
% for s = [1,2,6,7]
%     for ch =6:58
%         [~,apm(s,ch)] = ttest2(ac{s}(:,ch),ac{s}(:,59)); % media igual
%         [~,apd(s,ch)] = ttest(ac{s}(:,ch),ac{s}(:,59)); % media dif
%         [~,Dklpm(s,ch)] = ttest2(Dkl{s}(:,ch),Dkl{s}(:,59)); % media igual
%         [~,Dklpd(s,ch)] = ttest(Dkl{s}(:,ch),Dkl{s}(:,59)); % media dif
%         [~,j_pm(s,ch)] = ttest2(j_{s}(:,ch),j_{s}(:,59)); % media igual
%         [~,j_pd(s,ch)] = ttest(j_{s}(:,ch),j_{s}(:,59)); % media dif
%     end
%     figure
%     subplot(3,1,1); plot(apm(s,:)); hold on; plot(apd(s,:),'r'); title('p-values Acc');
%     subplot(3,1,2); plot(Dklpm(s,:)); hold on; plot(Dklpd(s,:),'r'); title('p-values Div');
%     subplot(3,1,3); plot(j_pm(s,:)); hold on; plot(j_pd(s,:),'r'); title('p-values Rayleight');
%
% end
%
% % type 2
% for s = [1,2,6,7]
%     for ch =6:58
%         [~,apm(s,ch)] = ttest(ac{s}(:,ch),ac{s}(:,59),'Alpha',0.01,'Tail','both'); % media igual
%         [~,apd(s,ch)] = ttest(ac{s}(:,ch),ac{s}(:,59),'Alpha',0.01,'Tail','right'); % media dif
%         [~,apd_(s,ch)] = ttest(ac{s}(:,ch),ac{s}(:,59),'Alpha',0.01,'Tail','left'); % media dif
%         [~,Dklpm(s,ch)] = ttest(Dkl{s}(:,ch),Dkl{s}(:,59),'Alpha',0.01,'Tail','both'); % media igual
%         [~,Dklpd(s,ch)] = ttest(Dkl{s}(:,ch),Dkl{s}(:,59),'Alpha',0.01,'Tail','right'); % media dif
%         [~,Dklpd_(s,ch)] = ttest(Dkl{s}(:,ch),Dkl{s}(:,59),'Alpha',0.01,'Tail','left'); % media dif
%         [~,j_pm(s,ch)] = ttest(j_{s}(:,ch),j_{s}(:,59),'Alpha',0.01,'Tail','both'); % media igual
%         [~,j_pd(s,ch)] = ttest(j_{s}(:,ch),j_{s}(:,59),'Alpha',0.01,'Tail','right'); % media dif
%         [~,j_pd_(s,ch)] = ttest(j_{s}(:,ch),j_{s}(:,59),'Alpha',0.01,'Tail','left'); % media dif
%     end
%     figure
%     subplot(3,1,1); plot(apm(s,:)); hold on; plot(apd(s,:),'r'); hold on; plot(apd_(s,:),'g'); title('p-values Acc');
%     subplot(3,1,2); plot(Dklpm(s,:)); hold on; plot(Dklpd(s,:),'r'); hold on; plot(Dklpd_(s,:),'g'); title('p-values Div');
%     subplot(3,1,3); plot(j_pm(s,:)); hold on; plot(j_pd(s,:),'r'); hold on; plot(j_pd_(s,:),'g'); title('p-values Rayleight');
%
% end
%
% % type 3
%
% for s = [1,2,6,7]
%     for ch =6:58
%         [hapm(s,ch),apm(s,ch)] = ttest(ac{s}(:,ch),ac{s}(:,59),'Alpha',0.01,'Tail','both'); % media igual
%         [hapd(s,ch),apd(s,ch)] = ttest(ac{s}(:,ch),ac{s}(:,59),'Alpha',0.01,'Tail','right'); % media dif
%         [hapd_(s,ch),apd_(s,ch)] = ttest(ac{s}(:,ch),ac{s}(:,59),'Alpha',0.01,'Tail','left'); % media dif
%         %         [~,Dklpm(s,ch)] = ttest(Dkl{s}(:,ch),Dkl{s}(:,59),'Alpha',0.01,'Tail','both'); % media igual
%         %         [~,Dklpd(s,ch)] = ttest(Dkl{s}(:,ch),Dkl{s}(:,59),'Alpha',0.01,'Tail','right'); % media dif
%         %         [~,Dklpd_(s,ch)] = ttest(Dkl{s}(:,ch),Dkl{s}(:,59),'Alpha',0.01,'Tail','left'); % media dif
%         %         [~,j_pm(s,ch)] = ttest(j_{s}(:,ch),j_{s}(:,59),'Alpha',0.01,'Tail','both'); % media igual
%         %         [~,j_pd(s,ch)] = ttest(j_{s}(:,ch),j_{s}(:,59),'Alpha',0.01,'Tail','right'); % media dif
%         %         [~,j_pd_(s,ch)] = ttest(j_{s}(:,ch),j_{s}(:,59),'Alpha',0.01,'Tail','left'); % media dif
%     end
%
%     figure
%     subplot(2,1,1); plot(6:58,mean(ac{s}(:,6:58)))%; hold on; plot(6:58,apd(s,6:end),'r'); hold on; plot(6:58,apd_(s,6:end),'g'); title('p-values Acc');
%     subplot(2,1,2); plot(6:58,apm(s,6:end)); hold on; plot(6:58,apd(s,6:end),'r'); hold on; plot(6:58,apd_(s,6:end),'g'); title('p-values Acc');
%     %xlim([6,58])
%     %     subplot(3,1,2); plot(6:58,Dklpm(s,6:end)); hold on; plot(6:58,Dklpd(s,6:end),'r'); hold on; plot(6:58,Dklpd_(s,6:end),'g'); title('p-values Div');
%     %xlim([6,58])
%     %     subplot(3,1,3); plot(6:58,j_pm(s,6:end)); hold on; plot(6:58,j_pd(s,6:end),'r'); hold on; plot(6:58,j_pd_(s,6:end),'g'); title('p-values Rayleight');
%     %xlim([6,58])
%
% end

% cc_ = [44,6,19,6]
% cc_ = I+5;
% ci= 1;
% for s = [1,2,6,7]
%     r_(ci,1) = mean(ac{s}(:,cc_(ci)));
%     r_(ci,2) = std(ac{s}(:,6));
%     ci = ci+1;
% end
% [mean(r_(:,1)),mean(r_(:,2))]
%
% hold on
% for s = [1,2,6,7]
%     plot(1:59,mean(ac{s}),'LineWidth',2)
% end
% legend('S1','S2','S6','S7')
%
% % plot
% figure
% posc =posch{1}(1:6);
% hold on
% scatter3(pos(:,1),pos(:,2),pos(:,3),'filled',jet(59));
% % scatter3(pos(posc,1),pos(posc,2),pos(posc,3),'MarkerEdgeColor',[0 0.1 1]);
% text(pos(:,1),pos(:,2),pos(:,3),electrodes);
%
% % for s = [3,4]
% % for i =1:4; j_{s}(i,:) = jray{s,i}; end
% % figure
% % errorbar(mean(j_{s}),std(j_{s}))
% % end
%
%
% xpru = (trace((W*Cov{2}*W')'*(W*Cov{1}*W')'))/(sqrt(trace(W*Cov{2}*W'*W*Cov{2}*W'))*sqrt(trace(W*Cov{1}*W'*W*Cov{1}*W')))
%
% norm(diag(W*Cov{1}*W')-diag(W*Cov{1}*W'))
%
% I=[39,1,14,1] % wang
% % I=[18,3,1,2]; % tp
% ch=I+5;
% set(0,'DefaultFigureWindowStyle','docked')
% figure(1)
% a=1;
% b=1;
% hold on
% for s = [1,2,6,7]
%     subplot(2,4,b)
%     % cl = gray(59);
%     % cl =flip(cl);
%     posc = posch{s}(1:ch(a));
%     % for ic =1:59
%     scatter3(pos(posc,1),pos(posc,2),pos(posc,3),300,'filled')%,'MarkerFaceColor',[cl(ic,:)]);
%     hold on
%     % end
%     view(0,90)
%     text(pos(:,1),pos(:,2),pos(:,3),electrodes);
%     a=a+1;
%     xlim([-1 1]); ylim([-1 1]);
%     b=b+1;
% end