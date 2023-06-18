
%% MAIN IWCSP
% COHEN2017: Using spatiotemporal source separation to identify
% prominent features in multichannel data without sinusoidal filters
%
% L.F. Velasquez-Martinez
% -------------------------------------------------------------------------
%% Initial parameters
clc, clear %all;            % parametros para manejo del codigo
% Chan = 8;  % canal (8-C3) - (10-Cz) - (12-C4) se puede colocar cualquiera de los 22.
tam = 500; % 280 = 1seg 500 = 2seg - ventana de analisis.
load('D:\Luisa\Codes\BCICIV_2a.mat');            % raw data
name = 'BCICIV2a';
ns = size(X,1);
sub = 1:ns;
Xd_ = cell(1,ns);
lb = cell(1,ns);
labels = [1];                          % clases
for i = 1:ns
    ind{i} = ismember(y{i},labels);
    Xd_{i} = X{i}(ind{i});               % belonging to labels
    lb{i} = y{i}(ind{i});
    Xdt{i} = X{i}(ind{i})';
end
fprintf('Loaded data: %s\n',name);


% %% Filter Band
% clc;
% a = [0.5 40];                                  % rango del filtro.
% % filtrado de cad sujeto.
% 
% for k = 1:9
%     fprintf(['Sujeto: ' num2str(k) ' de 9' '\n'])
%     Xfreq{k} = fcnfiltband(Xd_{k},fs,a,5);           % funcion Filtro.
%     %Xfreqt{k} = fcnfiltband(Xdt{k},fs,a,5);           % funcion Filtro.
% end
% Xd_ = Xfreq;
% %Xdt = Xfreqt;

%% Selecting data for analysis
tini = 2; % (MI segment)                             % tiempo de señal
if tam == 250; trfin = 4; else;  trfin =4; end
% samples = (trfin-tini)*fs;
Xddc = fncCutdata(Xd_,tini,trfin,fs);   % selecciona y filtra 1x500
% Xdc = cellfun(@(Xa) cellfun(@(X) X-mean(X),Xd_,'UniformOutput',false),Xdd,'UniformOutput',false);

trini = 0;  % (Reference segment)                      % tiempo de referencia.
if tam == 250;    trfin =2; else;  trfin = 2; end
Xdrc = fncCutdata(Xd_,trini,trfin,fs);                      % selecciona y filtra 1x500
% Xrc = cellfun(@(Xa) cellfun(@(X) X-mean(X),Xd_,'UniformOutput',false),Xdr,'UniformOutput',false);


%% Analysis 
for i = 1:ns
            
    Xtrd = Xddc{i};%Xdata  
    Xtrr = Xdrc{i};%Xreference
    
    Xdc1 = Xtrd(lb{i}==1);
    Xdc2 = Xtrd(lb{i}==2);
       
    % spatial pattern
%     Ws = fncCSPcohen(Xtrd,Xtrr,0,0); % Finding temporal filters
    Ws = fncCSPcohen(Xdc1,Xdc2,0,0); % Finding temporal filters
    [X] = fncrotate(Ws, Xdc1, Xdc2, Xtrr, Xdt{i}); %Rotating data
%     [X] = fncrotate(Ws, Xtrd, Xtrr, Xdt{i}); %Rotating data
%     Xdc = X{1}; Xrc = X{2}; Xdtr = X{3};
    Xdc1r_ = X{1}; Xdc2r_ = X{2}; Xrc = X{3}; Xdtr = X{4};
    Xdc(lb{i}==1) = Xdc1r_;  Xdc(lb{i}==2) = Xdc2r_; 
    
    %% hankelizar - time. 
    w = 350;
    Xdd = fncHankmod(Xdc,w); % data
    Xdr = fncHankmod(Xrc,w); %reference
    
    %% Etapa de diseño
    [Wt,eig] = fncCSPcohen(Xdd,Xdr,1,0); % Finding temporal filters
    figure(i)
%     [Wt,eig] = fncCSPcohen(Xdd,Xdr,1,0); % Finding temporal filters
%     subplot(1,3,1); plot(Wt'),title('Wt') 
%     subplot(1,3,2); fncFFT2(Wt(1,:),Wt(end,:),250); %ylim([10^(-4),1]);
%     xlim([0,60]);title('spectrum') 
% %     subplot(1,3,2); fncFFT(Wt(1,:),250); hold on; fncFFT(Wt(end,:),250); ylim([10^(-4),1]);xlim([0,60]);title('spectrum') 
%     subplot(1,3,3); stem(eig([1,end]));xlim([0.5,2.5]);title('eig') 
    figure(i)                    % plot Wt
    subplot(3,1,1)
    plot(Wt),title('Wt')         % plot fft de Wt.
    subplot(3,1,2)
    fncFFT(Wt,250);
    
    %% Yt = Wt*Xc
    Yt{i} = cellfun(@(X) conv(X,Wt,'same'),Xdtr','UniformOutput',false);
%     Yt{i} = cellfun(@(X) conv(X,Wt,'same'),Xdtr','UniformOutput',false);
    mYt{i} = mean(cell2mat(Yt{i}),1);
    sYt{i} = std(cell2mat(Yt{i}),1);

    
    %% Image  Mean - STD todos los sujetos
    set(0,'DefaultFigureWindowStyle','docked')

    t=0:1/250:7-1/250;
    MYT = mYt{i};     SYT = sYt{i};
    subplot(3,1,3)
    hold on
    plot(t',MYT,'b','LineWidth',2)
    plot([2 2],[min(MYT) max(MYT)],'--r','LineWidth',2)
    plot([3.25 3.25],[min(MYT) max(MYT)],'--g','LineWidth',2)
    plot([6 6],[min(MYT) max(MYT)],'--m','LineWidth',2)
    hold off
    title(['Sujeto: ' num2str(i)])

%     for tr = 1:20%numel(Yt{i})
%         figure; plot(0:1/250:7-1/250,cell2mat(Yt{i}(tr)))
%         pause(1);
%     end
%     pause();
%     close all
end
