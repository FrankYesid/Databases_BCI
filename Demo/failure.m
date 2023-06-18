% %%
clear all; clc; %close all;
% load('C:\Users\lfvelasquezm\Dropbox\ERD\Codes\TP\Matlab_con_csp\BCICIV_1.mat')   % carga base de datos EEG. BCICIV_1
% name  = 'BCIIV_1';
% fprintf('Loaded data: %s\n',name);
% ns = numel(X);
% sub = 1:ns;
% fs = 100;
% load('electrodesBCICIV1.mat')

load('C:\Users\lfvelasquezm\Dropbox\ERD\Codes\signal\BCICIV_2a.mat')   % carga base de datos EEG. BCICIV_1
ns = numel(X);
%parameters
tm = 2; % time embedding
w = tm*fs; 
trini = 0.0; % Inital reference time 
tini = 2.5; % Inital data time   

a = [4,36];

pos_ = [4,6,11,12,13,14,15,16,17,20,21,22,23,24,25,26,28,29,30,31,33,34,....
        35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,51,52,53,54,56,57,58,...
        59,60,61,62,65,66,67,68,69,70,71,76,78,85,87];

    
Wt = cell(1,ns);
freqWt = cell(1,ns);
Ch_img = cell(1,ns);

w_ = 250;


for ia = 1%: length(a)
    tic
    %filters
    aa = [a(ia),a(ia+1)+4];
    
    set(0,'DefaultFigureWindowStyle','docked')
    for s = [8]
%         tic
%         fprintf(['Sujeto: ' num2str(s) '\n'])
        Xc = X{s};
%         Xc = fcnfiltband(X{s},fs,aa,8);           % funcion Filtro.
        ciw = w/2+1;
        for iw = w_
             tic
%              Selecting class
%             labels = [1,2]                          % clases
%                  fprintf(['Class: ' num2str(labels) '\n'] )
%                  ind{s} = ismember(y{s},labels);
%                  Xc = X{s}(y{s}==labels);
%                  lb = y{s}(ind{s});
            
            
            for chan = 8%1: 59
                
                Xc1 = X{s}(y{s}==1);
                               
                Xd = cellfun(@(L) L(:,chan),Xc1,'UniformOutput',false);
                
                
                % Selecting data for analysis
                % (MI segment)                             % tiempo de señal
                tdfin = tini+tm;
                Xd1 = fcnfiltband(Xd,fs,[4,40],8);
                Xdc = fncCutdata_(Xd1,tini,tdfin,fs);   % selecciona y filtra 1x100
                for tr = 1: length(Xdc)
                    FXd(tr,:) = fncFFT(Xdc{tr},fs,0);
                end
                
                %clase 2
                Xc2 = X{s}(y{s}==2);
                Xd2_ = cellfun(@(L) L(:,chan),Xc2,'UniformOutput',false);
                % Selecting data for analysis
                Xd2_ = fcnfiltband(Xd2_,fs,[4,40],8);
                Xdc2 = fncCutdata_(Xd2_,tini,tdfin,fs);
                
                
                % (Reference segment)                      % tiempo de referencia.
                trfin = trini +tm;
                Xd2 = fcnfiltband(Xd,fs,[4,40],8);          
                Xrc = fncCutdata_(Xd2,trini,trfin,fs);                      % selecciona y filtra 1x500
                for tr = 1: length(Xrc)
                    FXr(tr,:) = fncFFT(Xrc{tr},fs,0);
                end
                
                %%
%                 T = 1/fs;           
%                 L = 500;             
%                 t = (0:L-1)*T;
%                 
%                 for i = 1 :100; Xdc{i} = 0.7*sin(2*pi*50*t)'; end
%                     
                % hankelizar - time. 1x500
                Xdd = fncHankeling_(Xdc,iw); % % data % fncHank_(Xdc,iw);%
                Xdr = fncHankeling_(Xrc,iw);%f %reference % fncHank_(Xrc,iw);%
                Xdd2 = fncHankeling_(Xdc2,iw);%f %reference % fncHank_(Xrc,iw);%
                
                [covXd1,~]= fncCov(Xdd,1,0);               % Data covariance
                [covXd2,~]= fncCov(Xdd2,1,0);               % Data covariance
                [covXr,~] = fncCov(Xdr,1,0);
                
                
                [Wtemp,eig] = fncCSPcohen(Xdd,Xdr,1,3); % Finding
                if sum(imag(Wtemp))~= 0
                    chimg(chan) = 1;
                else chimg(chan) = 0; 
                end
                
                Wtemp = Wtemp - mean(Wtemp);
                wtemp(chan,:) = Wtemp;
                %                 fprintf(['TP: Channel - ' num2str(chan) '\n'])
                %                 setmp(chan,:) = fncSamEnt(Wtemp,1,w,w/3);
                [freq(chan,:),~] = fncFFT(Wtemp,fs,0);
               
            end
            ciw = ciw-10;
            time = toc;
            fprintf('Embedding. %d. Time: %f\n',iw,time)
            Wt{s}{ia}{iw} = wtemp;
            wtemp = [];
            freqWt{s}{ia}{iw} = freq;
            freq = [];
            Ch_img{s}{ia}{iw} = chimg;
            %             SE{s}{ia}{labels} = setmp;
            
%             Xc = [];
            %             end
            
            %%
            
            %        plot
            
            %test
            
            %
            %         set(0,'DefaultFigureWindowStyle','docked')
            %         figure('Name',['Sub ' num2str(s) ' Temporal Patterns - Frequencies ' num2str(aa(1))...
            %             '-' num2str(aa(2))],'NumberTitle','off')
            %         for chan = 1: 59
            %             subplot(10,9,pos_(chan))
            %             plot(0:1/fs:tm-1/fs,Wt{ia}{s}(chan,:))
            %             if chan == 27 || chan == 29 || chan == 31
            %                 plot(0:1/fs:tm-1/fs,Wt{ia}{s}(chan,:),'r')
            %             end
            %
            %             %title(num2str(electrodes{chan}),'FontSize',16,'Interpreter','latex')
            %             ylim([-0.2,0.2]);
            %         end
            %
            %         figure('Name',['Sub ' num2str(s) ' Spectral - Frequencies ' num2str(aa(1))...
            %             '-' num2str(aa(2))],'NumberTitle','off')
            %         f = fs*(0:(w/2))/w;
            %         for chan = 1: 59
            %             subplot(10,9,pos_(chan))
            %             freq=fncFFT(Wt{ia}{s}(chan,:),fs);
            %             freq = freq./max(freq);
            %             plot(f,freq)
            %             if chan == 27 || chan == 29 || chan == 31
            %                 plot(f,freq,'r')
            %             end
            %             xlim([0,40])
            %             %title(num2str(electrodes{chan}),'FontSize',16,'Interpreter','latex')
            %         end
            
            
        end
%         time = toc;
%         fprintf('Sub. %d. Freq. %d. Class. %d Time: %f\n',s,ia,labels,time)
        
    end
end


for s = [7]
    figure
    for ia = 1
        for iw = w_
            f = fs*(0:(iw/2))/iw;
            fm = freqWt{s}{1}{iw};
            plot(f,mean(fm));
            hold on
        end
        legend('w=100','w=110','w=120','w=130','w=140','w=150','w=160','w=170','w=180','w=190')
    end
end

% f = fs*(0:(w/2))/w;
% for s = [2,7]
%     for cl = [1,2]
% %         figure('Name',['Sub ' num2str(s) ' Spectral - Clase ' num2str(cl)],'NumberTitle','off')
%         for ia = 1%:length(a)
%             aa = [a(ia),a(ia)+4];
%             fm = freqWt{s}{ia}{cl};
%             errorbar(f,mean(fm),std(fm));
%             %l{ia}= {([num2str((ia)) '-' num2str((ia)+4)])};
%             hold on
%         end
%         legend({'4-30'})
% %         legend({'4-8','6-10','8-12','10-14','12-16','14-18','16-20','18-22',...
% %             '20-24','22-26','24-28','26-30','28-32','30-34','32-38','34-38','36-40'})
%         
%     end
% end