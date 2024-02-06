clear all;
close all;
%% ambos labels
% inda = randperm(60);
inda = 1:60;
%%
dB = [-10];
sub = 1;
fs = 100;
X = cell(sub,1);
y = cell(sub,1);
fre1 = 8;

t = 2:1/fs:6-1/fs;
indx = randperm(25);
tr_ = ones(1,60);

for s = 1
    a = 1; X{s} = cell(1,60);
    for tr = inda
 
        if a <=30
            cla = 1;
        else
            cla = 2;
        end
        for ch = 1
            %                             t1 = wgn(200,1,dB)';
            %                             t3 = wgn(200,1,dB)';
            %                             if cla == 1
            %                                 A1 = 20;
            %                                 A2 = round(A1/10*rand);
            %                                 t2 = A1*sin(2*pi*fre1*t) + A2*sin(2*pi*fre2*t);
            %                             elseif cla ==  2
            %                                 A2 = 20;
            %                                 A1 = round(A2/10*rand);
            %                                 t2  = A1*sin(2*pi*fre1*t) + A2*sin(2*pi*fre2*t);
            %                             end
            %                             X{s}{1,tr}(:,ch) = [t1,t2,t3];
            %                             y{s}(1,tr) = cla;            
            tf = wgn(800,1,dB)';
            if cla == 1
                %                 A1 = 20;
                %                 A2 = round(A1/10*rand);
                t2 = 20*sin(2*pi*fre1*t);%20*chirp(t,8,400/fs,24);% %20*sin(2*pi*fre1*t)+
%                 t2 = t2(end:-1:1);
%                 figure; subplot(2,1,1); plot(t2);
%                 subplot(2,1,2); fncFFT(t2,fs);
                %t2 = A1*sin(2*pi*fre1*t) + A2*sin(2*pi*fre2*t);
                tt = zeros(1,200);
                sig{1} = tt(1:end-indx(tr_(tr)));
                sig{2} = t2;
                sig{3} = tt(end-indx(tr_(tr))+1:end);
                sig{4} = tt;
                sig_ = cell2mat(sig);
                X{s}{1,tr}(:,ch) = (tf)+sig_;%tf+[tt,t2,tt];
                y{s}(1,tr) = cla;
            elseif cla ==  2
                %                 A2 = 20;
                %                 A1 = round(A2/10*rand);
                %                 t2  = A1*sin(2*pi*fre1*t) + A2*sin(2*pi*fre2*t);
                t2 = zeros(1,length(t));
                tt = zeros(1,200);
                X{s}{1,tr}(:,ch) = tf+[tt,t2,tt];
                y{s}(1,tr) = cla;
            end
            
        end
        a = a+1;
    end
end


% X{1} = fcnfiltband(X{1},fs,[8,12],15);   

% TP
ns = numel(X);
sub = 1:ns;

fs = 100;
Xddd = cell(1,ns);
lb = cell(1,ns);
N_channel = size(X{1}{1},2);
labels = [1,2];
for i = 1:ns
    ind{i} = ismember(y{i},labels);
    Xddd{i} = X{i}(ind{i}); %belonging to labels
    lb{i} = y{i}(ind{i});
end

for ch = 1:N_channel
    tic
    Chan =ch;
    Xa_ = cell(1,ns);                             % celda almacenadora de solo los canales.
    for i = 1:ns
        Xa_{i} = cellfun(@(L) L(:,Chan),Xddd{i},'UniformOutput',false);
    end
    % Xa_ = cellfun(@(Xa) cellfun(@(X) X-mean(X),Xa,'UniformOutput',false),Xa_,'UniformOutput',false);
    % Xd_ = Xa_; clear Xa_;
    
    % Selecting data for analysis
    tini = 2; % (MI segment)                             % tiempo de señal
    trfin = 3.5;
    
    samples = (trfin-tini)*fs;
    Xdd = fncCutdata(Xa_,tini,trfin,fs);   % selecciona y filtra 1x100
    %     Xdc = cellfun(@(Xa) cellfun(@(X) X-mean(X),Xa_,'UniformOutput',false),Xdd,'UniformOutput',false);
    Xdc = Xdd;
    
    trini = 0;  % (Reference segment)                      % tiempo de referencia.
    trfin = 1.5;
    samplesr = (trfin-trini)*fs;
    Xdr = fncCutdata(Xa_,trini,trfin,fs);                      % selecciona y filtra 1x500
    %     Xrc = cellfun(@(Xa) cellfun(@(X) X-mean(X),Xa_,'UniformOutput',false),Xdr,'UniformOutput',false);
    Xrc = Xdr;
    
    % hankelizar - time. 1x500
    w = 1.5*fs;
    Xdd_ = fncHank(Xdc,w); % data
    Xdr_ = fncHank(Xrc,w); %reference
    
    % Analysis
    for i = 1 %[1 2 3 4 5 6 7 8 9]
        
        %
        % Regularization
        %Xtrd = cellfun(@(X) X+(0.1.*eye(250,500).*X),Xtrd,'UniformOutput',false);
        %Xtrr = cellfun(@(X) X+(0.1.*eye(250,500).*X),Xtrr,'UniformOutput',false);
        
        % Etapa de diseño
        Wt = fncCSPcohen(Xdd_{1},Xdr_{1},2,3); % Finding temporal filters $%(Xd,Xdr,mode,flag)
        %         wwtt{i}{ch} = Wt;
        %         figure(i)                    % plot Wt
        %         subplot(3,1,1)
        %         plot(Wt),title('Wt')         % plot fft de Wt.
        %         subplot(3,1,2)
        %         PP1{i}{ch} = fncFFT(Wt,250);
        
        %         %% Yt = Wt*Xc
        %         Yt_conv = cellfun(@(X) conv(X,Wt','same'),Xa_{i},'UniformOutput',false);
        %         Yt{i}{ch} = cellfun(@(X) envelope(X)',Yt_conv,'UniformOutput',false);
        %     save('Yt','Yt')              % guardar matriz
        %     figure; plot(Yt{i}(1))
    end
    time = toc;
    fprintf(['Canal: ' num2str(ch) ' Time: ' num2str(time) '\n'])
end

% plot
set(0,'DefaultFigureWindowStyle','docked')

figure
subplot(3,1,1)
plot(0:1/fs:8-1/fs,X{1}{1},'LineWidth',1.5)
hold on 
plot(0:1/fs:8-1/fs,X{1}{31},'r','LineWidth',1.5)
xlim([1 3])
title(['Signal - Class ' num2str(labels) ' in ' num2str(dB) ' dB '])
xlabel('Time (s)','Interpreter','latex')
ylim([-20 20])
subplot(3,1,2)
plot(0:1/fs:2-1/fs,Wt)
xlim([1,3]);

title('Temporal Pattern','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex')
subplot(3,1,3)
f = fs*(0:(w/2))/w;
freq=fncFFT(Wt,fs);
plot(f,freq)
% spectrogram(Wt,'yaxis')
% [s,f,t,p] = spectrogram(Wt,[],[],[],100);
%spectrogram(Wt,50,25,[],100,'yaxis');
%colorbar('off')
% xlim([1,3])
% ylim([0,50])
% title('Frecuency','Interpreter','latex')
% xlabel('Frequency (Hz)','Interpreter','latex')

% figure
% for i =1 :60
%     plot(0:1/fs:8-1/fs,X{1}{i}); hold on
% end
% plot()