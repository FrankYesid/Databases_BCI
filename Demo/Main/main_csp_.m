% size%% MAIN CSP
clc, clear %all;*
%% parameters
% dt = 1,2,3:  1:BCICIV_2a
%%
fs =100;
% for dt = 1:3

% dt = 1; % load data
% dataset = dt;
% [Xdr,lb,cv,sub,name,ind] = fncloaddata(dataset);
load('F:\Database\BCICIV_1\dataBase.mat');          % raw data
load('Datos.mat')                                   % data pos.
name = 'BCICIV_1';
ns = size(X,1);
sub = 1:ns;
Xd_ = cell(1,ns);
lb = cell(1,ns);
labels = [1 2];
for i = 1:ns
    ind{i} = ismember(y{i},labels);
    Xd_{i} = X{i}(ind{i}); %belonging to labels
    lb{i} = y{i}(ind{i});
end

%% (MI segment)
tini = 2.5;
tfin = 4.5;
Xdr = fncCutdataf(Xd_,tini,tfin,fs,[2,40]);  % selecciona y filtra
fprintf('Loaded data: %s\n',name);

%% Analysis
np = 5; % Kfolds
[S,F] = ndgrid(sub,1:np);
avtrain = zeros(size(S));
avtest = zeros(size(S));

for i = 1:numel(S)
    tic
    s = S(i);
    f = F(i);
    
    Xd = Xdr{s};
    
    tr_ind = cv{s}.training(f); tr_ind = tr_ind(ind{s});
    ts_ind = cv{s}.test(f); ts_ind = ts_ind(ind{s});
    
    Xtrain = Xd(tr_ind);
    Xtest = Xd(ts_ind);
    idtr = find(tr_ind==1);
    ltr = lb{s}(idtr)';
    idts = find(ts_ind==1);
    lts = lb{s}(idts)';
    
    %     for ch = su
    %     for chan = 6:59
    %% CSP
    W = fncCSP(Xtrain,ltr,3);
    %A = inv(W);
    %feature
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
    tiempo = toc;
    fprintf('Exp - %d. Time: %f\n',i,tiempo)
    fprintf('Acc    : Training %f. Test %f\n', avtrain(i)*100,avtest(i)*100)
    %[i,k]
    %     end
end

%% CSP
% fprintf('Testing Results: %s \n', name)
% for i = 1:numel(sub)
%     acc_s(i,1) = mean(avtest(i,:,:))*100;
%     acc_s(i,2) = std(avtest(i,:,:))*100;
%     fprintf('Acc : %f. +- %f\n', acc_s(i,1),acc_s(i,2))
% end
%%
% save(strcat('Results\',sprintf('ResultsCSP_%s.mat',name)))
% end

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
