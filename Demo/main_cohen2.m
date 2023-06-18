
%% MAIN IWCSP
% COHEN2017: Using spatiotemporal source separation to identify
% prominent features in multichannel data without sinusoidal filters
%
% L.F. Velasquez-Martinez
% -------------------------------------------------------------------------

clc, clear %all;            % parametros para manejo del codigo
Chan = 8;  % canal (8-C3) - (10-Cz) - (12-C4) se puede colocar cualquiera de los 22.
YYtt = 1;  % 0 - Yt = Wt.Xc, 1 - Yt = Wt*Xc
tam = 500; % 280 = 1seg 500 = 2seg - ventana de analisis.

%%
load('D:\Archivos\Databasses_BCI\BCICIV_2a.mat');            % raw data
name = 'BCICIV2a';
ns = size(X,1);
sub = 1:ns;
Xd_ = cell(1,ns);
lb = cell(1,ns);
labels = [2];                          % clases
for i = 1:ns
    ind{i} = ismember(y{i},labels);
    Xd_{i} = X{i}(ind{i});               % belonging to labels
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
        Xfreq{k} = fcnfiltband(X{k},fs,a,5);           % funcion Filtro.
    end
    X = Xfreq;
end

for ch = 1:22
    %% Channels C3, C4
    tic
    Chan =ch;
    Xa_ = cell(1,ns);                             % celda almacenadora de solo los canales.
    
    % load revchan.mat                            % relevancia Renyi entropia
    % for i =1:ns
    %     Chan(i) = find(mat(i,:) == max(mat(i,:)));
    % end
    %
    % load relevance.mat                          % relevancia.
    % class = 1;
    % for i = 1:ns
    %     Chan(i) = find(relevance.R{i,class} == max(relevance.R{i,class}));
    % end
    
    for i = 1:ns
        Xa_{i} = cellfun(@(L) L(:,Chan),X{i},'UniformOutput',false);
    end
    Xa_ = cellfun(@(Xa) cellfun(@(X) X-mean(X),Xa,'UniformOutput',false),Xa_,'UniformOutput',false);
    % Xd_ = Xa_; clear Xa_;
    
    %% Selecting data for analysis
    tini = 2; % (MI segment)                             % tiempo de señal
    if tam == 250
        trfin = 4;
    else
        trfin =4;
    end
    samples = (trfin-tini)*fs;
    Xdd = fncCutdata(Xa_,tini,trfin,fs);   % selecciona y filtra 1x500
    % for s = 1: ns
    %     Xdc{s} = cellfun(@(X) X-mean(X),Xdd{s},'UniformOutput',false);
    % end
    Xdc = cellfun(@(Xa) cellfun(@(X) X-mean(X),Xa,'UniformOutput',false),Xdd,'UniformOutput',false);
    
    trini = 0;  % (Reference segment)                      % tiempo de referencia.
    if tam == 250
        trfin =2;
    else
        trfin = 2;
    end
    samplesr = (trfin-trini)*fs;
    Xdr = fncCutdata(Xa_,trini,trfin,fs);                      % selecciona y filtra 1x500
    % for s = 1: ns
    %     Xrc{s} = cellfun(@(X) X-mean(X),Xdr{s},'UniformOutput',false);
    % end
    Xrc = cellfun(@(Xa) cellfun(@(X) X-mean(X),Xa,'UniformOutput',false),Xdr,'UniformOutput',false);
    
    %% hankelizar - time. 1x500
    w = 350;
    Xdd = fncHank(Xdc,w); % data
    Xdr = fncHank(Xrc,w); %reference
    
    %% Analysis
    for i = [1 2 3 4 5 6 7 8 9]
        
        Xtrd = Xdd{i};%Xdata  ----
        Xtrr = Xdr{i};%Xreference
        
        %% Regularization
        %Xtrd = cellfun(@(X) X+(0.1.*eye(250,500).*X),Xtrd,'UniformOutput',false);
        %Xtrr = cellfun(@(X) X+(0.1.*eye(250,500).*X),Xtrr,'UniformOutput',false);
        
        %% Etapa de diseño
        Wt = fncCSPcohen(Xtrd,Xtrr,1,0); % Finding temporal filters
        %         figure(i)                    % plot Wt
        %         subplot(3,1,1)
        %         plot(Wt),title('Wt')         % plot fft de Wt.
        %         subplot(3,1,2)
        %         fncFFT(Wt,250);
        
        %% Yt = Wt*Xc
        Yt_conv = cellfun(@(X) conv(X,Wt','same'),Xa_{i},'UniformOutput',false);
        Yt{i}{ch} = cellfun(@(X) envelope(X)',Yt_conv,'UniformOutput',false);
        %     save('Yt','Yt')              % guardar matriz
        %     figure; plot(Yt{i}(1))
        
        %% Mean - STD
        %     load('Yt.mat')
        %     for j = 1:numel(Yt{i})
        
        %         mYt{i} = mean(cell2mat(Yt{i}),1);
        %         sYt{i} = std(cell2mat(Yt{i}),1);
        %     end
        %     mYt = cellfun(@(Xa) mean(cell2mat(Xa),1),Yt,'UniformOutput',false);
        %     sYt = cellfun(@(Xa) std(cell2mat(Xa),1),Yt,'UniformOutput',false);
        %     save('mYt.mat','mYt')        % guardar matriz de mean - std
        %     save('sYt.mat','sYt')
        
        %% Image  Mean - STD todos los sujetos
        %         set(0,'DefaultFigureWindowStyle','docked')
        %         if Chan == 8
        %             nam = 'C3';
        %         else
        %             nam = 'C4';
        %         end
        %     load('mYt.mat')
        %     load('sYt.mat')
        %         t=0:1/250:7-1/250;
        %         MYT = mYt{i};
        %         SYT = sYt{i};
        %         subplot(3,1,3)
        %         hold on
        %         plot(t',MYT,'b','LineWidth',2)
        %         plot([2 2],[min(MYT) max(MYT)],'--r','LineWidth',2)
        %         plot([3.25 3.25],[min(MYT) max(MYT)],'--g','LineWidth',2)
        %         plot([6 6],[min(MYT) max(MYT)],'--m','LineWidth',2)
        %         hold off
        %         title(['Sujeto: ' num2str(i) ' Canal: ' nam])
        
        %     for tr = 1:20%numel(Yt{i})
        %         figure; plot(0:1/250:7-1/250,cell2mat(Yt{i}(tr)))
        %         pause(1);
        %     end
        %     pause();
        %     close all
    end
    time = toc;
    fprintf(['Canal: ' num2str(ch) ' Time: ' num2str(time) '\n'])
end

nat = 0; %% para ver solo las areas
if nat == 1
    x = 1 : 1750;
    curve1 = a';
    curve2 = b';
    plot(x, curve1, 'c', 'LineWidth', 2);alpha(0.25)
    hold on;
    plot(x, curve2, 'c', 'LineWidth', 2);alpha(0.25)
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    fill(x2, inBetween, 'c');
    alpha(0.25)
end
%%
load Yt1.mat
load Yt2.mat
clc;
% Filter
fil = 0;
if fil == 1
    for as=1:2
        if as == 1
            X = Yt1;
        else
            X = Yt2;
        end
        % Filter Band
%         clc;
        a = [13 30];                                  % rango del filtro.
%         A = cell(9,1);                               % celdas de la señal filtrada.
        filt = 0;
        % filtrado de cad sujeto.
        for k = 1:9
            fprintf(['Sujeto: ' num2str(k) ' de 9' '\n'])
            Xfreq{k} = fcnfiltband2(X{k},fs,a,5);           % funcion Filtro.
        end
        if as == 1
            Yt1 = Xfreq;
        else
            Yt2 = Xfreq;
        end
        
    end
end


grap = 2;   % tipo de grafica 2- 2D o 3- 3D
if grap == 2    %%%%%%%%%%%% Headplot for channel ciclo spectrogram - 2D %%%%%%%%%%% %%
    set(0,'DefaultFigureWindowStyle','docked')                                  % grafica cada figura en una misma pantalla.
    posi = [4 9 10 11 12 13 15 16 17 18 19 20 21 23 24 25 26 27 31 32 33 39];   % ubicación de cada una de las graficas.   % posición de las graficas
    % rango de amplitud
    %     Lrang = [-0.8 0.8];                                                   % rango de colocar la señal.
    for sub = 3                                                              % graficador en 2D todos los canales
        i = 1;
        for ch = 1:22
            figure(sub)
            subplot(6,7,posi(i))
            xpru = [cell2mat(Yt1{sub}{ch});cell2mat(Yt2{sub}{ch})];
            m = max(xpru');
            imagesc(0:1/250.0:7-1/250,1:numel(m),...
                squeeze(bsxfun(@times,xpru,(1./m)')))                     %%
            if ch == 1:22
                title('C3')
            elseif ch == 12
                title('C4')
            else
                title([num2str(ch)])
            end
            %             xlim([1 t(end)])
            %             ylim([1 40]);
            yticks([])
            xticks([])
            if ch == 22
                xlabel('Time(seg)','Interpreter','latex')
                ylabel('Trials','Interpreter','latex')
                xticks([0 1 2 3 4 5 6])
            end
            hold on
            plot([2 2],[1 size(m,2)],'--r','LineWidth',2)
            hold off
            colormap jet
            drawnow
            i = i+1;
        end
        
    end
elseif grap ==3    %%%%%%%%%%%% Headplot for channel ciclo spectrogram - 3D %%%%%%%%%%% %%%
    set(0,'DefaultFigureWindowStyle','docked')                                  % grafica cada figura en una misma pantalla.
    posi = [4 9 10 11 12 13 15 16 17 18 19 20 21 23 24 25 26 27 31 32 33 39];   % ubicación de cada una de las graficas.    % posición de las graficas
    for sub = 3                                                              % graficador en 3D todos los canales
        i = 1;
        
        for ch = 1
            figure(sub+20)
            subplot(6,7,posi(i))
            surf(0:1/250.0:7-1/250.0,1:numel(Yt{sub}{ch}),squeeze(cell2mat(Yt{sub}{ch})))
            %         title([labels{ch}])
            if ch == 8
                title('C3')
            elseif ch == 12
                title('C4')
            else
                title([num2str(ch)])
            end
            xlim([0 7])
            ylim([0 numel(Yt{sub}{ch})])
            zticks([])
            shading interp
            view(-18,38)
            colormap jet
            grid on
            %                 hold on
            %                 %             plot3(t.*0+2,t.*5,t.*0,'-r','LineWidth',1.5)
            %                 %             plot3(t.*0+4.5,t.*5,t.*0,'-r','LineWidth',1.5)
            %                 hold off
            drawnow
            i = i+1;
        end
        
    end
end