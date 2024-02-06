clear
clc
%%
inda = randperm(60);
%%
dB = [20];
sub = 1;
fs = 100;
X = cell(sub,1);
y = cell(sub,1);
fre1 = 8;
fre2 = 22;
t = 2:1/fs:6-1/fs;

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
                A1 = 20;
                A2 = round(A1/10*rand);
                t2 = A1*sin(2*pi*fre1*t) + A2*sin(2*pi*fre2*t);
            elseif cla ==  2
                A2 = 20;
                A1 = round(A2/10*rand);
                t2  = A1*sin(2*pi*fre1*t) + A2*sin(2*pi*fre2*t);
            end
            tt = zeros(1,200);
            X{s}{1,tr}(:,ch) = tf+[tt,t2,tt];
            y{s}(1,tr) = cla;
        end
        a = a+1;
    end
end

%% TP
ns = numel(X);
sub = 1:ns;

fs = 100;
Xddd = cell(1,ns);
lb = cell(1,ns);
N_channel = size(X{1}{1},2);
for labels = 1:2
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
    trfin =4;
    
    samples = (trfin-tini)*fs;
    Xdd = fncCutdata(Xa_,tini,trfin,fs);   % selecciona y filtra 1x100
    %     Xdc = cellfun(@(Xa) cellfun(@(X) X-mean(X),Xa_,'UniformOutput',false),Xdd,'UniformOutput',false);
    Xdc = Xdd;
    
    trini = 0;  % (Reference segment)                      % tiempo de referencia.
    trfin = 2;
    samplesr = (trfin-trini)*fs;
    Xdr = fncCutdata(Xa_,trini,trfin,fs);                      % selecciona y filtra 1x500
    %     Xrc = cellfun(@(Xa) cellfun(@(X) X-mean(X),Xa_,'UniformOutput',false),Xdr,'UniformOutput',false);
    Xrc = Xdr;
    
    % hankelizar - time. 1x500
    w = 200;
    Xdd = fncHank(Xdc,w); % data
    Xdr = fncHank(Xrc,w); %reference
    
    % Analysis
    for i = 1 %[1 2 3 4 5 6 7 8 9]
        
        %
        % Regularization
        %Xtrd = cellfun(@(X) X+(0.1.*eye(250,500).*X),Xtrd,'UniformOutput',false);
        %Xtrr = cellfun(@(X) X+(0.1.*eye(250,500).*X),Xtrr,'UniformOutput',false);
        
        % Etapa de diseño
        Wt = fncCSPcohen(Xdd{i},Xdr{i},2,0); % Finding temporal filters
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

if labels == 1
    dd = [1,3,5];
else
    dd = [2,4,6];
end
subplot(3,2,dd(1))
plot(0:1/fs:8-1/fs,X{1}{1},'LineWidth',1.5)
title(['Signal - Class ' num2str(labels) ' in ' num2str(dB) ' dB '])
xlabel('Time (s)','Interpreter','latex')
ylim([-20 20])
subplot(3,2,dd(2))
plot(0:1/fs:2-1/fs,Wt)
title('Temporal Pattern','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex')
subplot(3,2,dd(3))
L = 200; Fs = 100;
f = Fs*(0:(L/2))/L;
freq=fncFFT(Wt,fs);
plot(f,freq./max(freq))
title('Frecuency','Interpreter','latex')
xlabel('Frequency (Hz)','Interpreter','latex')
end

%% ambos labels
inda = randperm(60);
%%
dB = [20];
sub = 1;
fs = 100;
X = cell(sub,1);
y = cell(sub,1);
fre1 = 8;
fre2 = 22;
t = 2:1/fs:6-1/fs;

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
                A1 = 20;
                A2 = round(A1/10*rand);
                t2 = A1*sin(2*pi*fre1*t) + A2*sin(2*pi*fre2*t);
            elseif cla ==  2
                A2 = 20;
                A1 = round(A2/10*rand);
                t2  = A1*sin(2*pi*fre1*t) + A2*sin(2*pi*fre2*t);
            end
            tt = zeros(1,200);
            X{s}{1,tr}(:,ch) = tf+[tt,t2,tt];
            y{s}(1,tr) = cla;
        end
        a = a+1;
    end
end

%% TP
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
    trfin =4;
    
    samples = (trfin-tini)*fs;
    Xdd = fncCutdata(Xa_,tini,trfin,fs);   % selecciona y filtra 1x100
    %     Xdc = cellfun(@(Xa) cellfun(@(X) X-mean(X),Xa_,'UniformOutput',false),Xdd,'UniformOutput',false);
    Xdc = Xdd;
    
    trini = 0;  % (Reference segment)                      % tiempo de referencia.
    trfin = 2;
    samplesr = (trfin-trini)*fs;
    Xdr = fncCutdata(Xa_,trini,trfin,fs);                      % selecciona y filtra 1x500
    %     Xrc = cellfun(@(Xa) cellfun(@(X) X-mean(X),Xa_,'UniformOutput',false),Xdr,'UniformOutput',false);
    Xrc = Xdr;
    
    % hankelizar - time. 1x500
    w = 200;
    Xdd = fncHank(Xdc,w); % data
    Xdr = fncHank(Xrc,w); %reference
    
    % Analysis
    for i = 1 %[1 2 3 4 5 6 7 8 9]
        
        %
        % Regularization
        %Xtrd = cellfun(@(X) X+(0.1.*eye(250,500).*X),Xtrd,'UniformOutput',false);
        %Xtrr = cellfun(@(X) X+(0.1.*eye(250,500).*X),Xtrr,'UniformOutput',false);
        
        % Etapa de diseño
        Wt = fncCSPcohen(Xdd{i},Xdr{i},2,0); % Finding temporal filters
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

subplot(3,1,1)
plot(0:1/fs:8-1/fs,X{1}{1},'LineWidth',1.5)
title(['Signal - Class ' num2str(labels) ' in ' num2str(dB) ' dB '])
xlabel('Time (s)','Interpreter','latex')
ylim([-20 20])
subplot(3,1,2)
plot(0:1/fs:2-1/fs,Wt)
title('Temporal Pattern','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex')
subplot(3,1,3)
L = 200; Fs = 100;
f = Fs*(0:(L/2))/L;
freq=fncFFT(Wt,fs);
plot(f,freq./max(freq))
title('Frecuency','Interpreter','latex')
xlabel('Frequency (Hz)','Interpreter','latex')

