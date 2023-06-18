%%
clear all; clc; %close all;
load('C:\Users\lfvelasquezm\Dropbox\ERD\Codes\TP\Matlab_con_csp\BCICIV_1.mat')   % carga base de datos EEG. BCICIV_1
name  = 'BCIIV_1';
fprintf('Loaded data: %s\n',name);
ns = numel(X);
sub = 1:ns;
fs = 100;
load('electrodesBCICIV1.mat')
%parameters
tm = 1.5; % time embedding
w = tm*fs; 
trini = 0.0; % Inital reference time 
tini = 2; % Inital data time   

a = 4:2:26;

pos_ = [4,6,11,12,13,14,15,16,17,20,21,22,23,24,25,26,28,29,30,31,33,34,....
        35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,51,52,53,54,56,57,58,...
        59,60,61,62,65,66,67,68,69,70,71,76,78,85,87];

for ia = 1: length(a)
    %filters
    aa = [a(ia),a(ia)+4];

    set(0,'DefaultFigureWindowStyle','docked')
    for s = 7%[1,2,6,7]
        tic
        fprintf(['Sujeto: ' num2str(s) ' de 9' '\n'])
        X{s} = fcnfiltband(X{s},fs,aa,8);           % funcion Filtro.
        
        for chan = 1: 59
            
            Xd = cellfun(@(L) L(:,chan),X{s},'UniformOutput',false);
            
            % Selecting data for analysis
            % (MI segment)                             % tiempo de señal
            trfin = tini+tm;
            Xdc = fncCutdata_(Xd,tini,trfin,fs);   % selecciona y filtra 1x100
            
            % (Reference segment)                      % tiempo de referencia.
            trfin = trini +tm;
            Xrc = fncCutdata_(Xd,trini,trfin,fs);                      % selecciona y filtra 1x500
            
            % hankelizar - time. 1x500
            Xdd = fncHank_(Xdc,w); % data
            Xdr = fncHank_(Xrc,w); %reference
            Wtemp = fncCSPcohen(Xdd,Xrc,2,3); % Finding
            Wtemp = Wtemp - mean(Wtemp);
            Wt{ia}{s}(chan,:) = Wtemp;
            
        end
%     end
% end
%    
% save('BCICIV1_TP_freqs.mat','Wt')

% for ia = 2: length(a)
%     %filters
%     aa = [a(ia),a(ia)+4];
%     
%     for s = 7%[1,2,6,7]
        set(0,'DefaultFigureWindowStyle','docked')
        figure('Name',['Sub ' num2str(s) ' Temporal Patterns - Frequencies ' num2str(aa(1))...
            '-' num2str(aa(2))],'NumberTitle','off')
        for chan = 1: 59
            subplot(10,9,pos_(chan))
            plot(0:1/fs:tm-1/fs,Wt{ia}{s}(chan,:))
            if chan == 27 || chan == 29 || chan == 31
                plot(0:1/fs:tm-1/fs,Wt{ia}{s}(chan,:),'r')
            end
            
            %title(num2str(electrodes{chan}),'FontSize',16,'Interpreter','latex')
            ylim([-0.2,0.2]);
        end
        
        figure('Name',['Sub ' num2str(s) ' Spectral - Frequencies ' num2str(aa(1))...
            '-' num2str(aa(2))],'NumberTitle','off')
        f = fs*(0:(w/2))/w;
        for chan = 1: 59
            subplot(10,9,pos_(chan))
            freq=fncFFT(Wt{ia}{s}(chan,:),fs);
            freq = freq./max(freq);
            plot(f,freq)
            if chan == 27 || chan == 29 || chan == 31
                plot(f,freq,'r')
            end
            xlim([0,40])
            %title(num2str(electrodes{chan}),'FontSize',16,'Interpreter','latex')
        end
        
        %         figure
        %         subplot(3,1,1)
        %         Xd1 = Xd(y{s}==1);
        %         Xd2 = Xd(y{s}==2);
        %
        %         for ii = 1:100
        %             %         plot(0:1/fs:8-1/fs,Xd1{ii},'LineWidth',1.5)
        %             plot(0:1/fs:8-1/fs,Xd1{ii},'LineWidth',1.5)
        %             hold on
        %         end
        %
        %         for ii = 1:100
        %             %         plot(0:1/fs:8-1/fs,Xd2{ii},'LineWidth',1.5)
        %             plot(0:1/fs:8-1/fs,Xd2{ii},'LineWidth',1.5)
        %             hold on
        %         end
        %         title(['Trials Sub ' num2str(s) ' - Class 1-2'],'Interpreter','latex')
        %         xlabel('Time (s)','Interpreter','latex')
        %         %     ylim([-20 20])
        %
        %         subplot(3,1,2)
        %         %     plot(0:1/fs:2-1/fs,Wt{s})
        %         plot(0:1/fs:tm-1/fs,Wt{s})
        %         xlim([0,8]);
        %         title('Temporal Pattern','Interpreter','latex')
        %         xlabel('Time (s)','Interpreter','latex')
        %
        %         subplot(3,1,3)
        %         %     f = fs*(0:(w/2))/w;
        %         %freq=fncFFT(Wt,fs);
        %         %     plot(f,freq./max(freq))
        %         %     spectrogram(Wt{s},'yaxis')
        %         % [s,f,t,p] = spectrogram(Wt,[],[],[],100);
        %         spectrogram(Wt{s},50,25,[],100,'yaxis');
        %         colorbar('off')
        %         xlim([0,8]);
%         time = toc;
%         fprintf('Suj: %d. Time: %f\n',s,time)
    end
    
end



