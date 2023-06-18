%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc
load('F:\Database\BCICIV_2a\BCICIV_2a.mat')
% load('G:\Dropbox\ERD\Codes\TP\Matlab_con_csp\BCICIV_1.mat')
% fs = 100;
Freq = [4,40];
subs = 8;
for s = subs
    temp{s} = fcnfiltband(X{s},fs,Freq,8);
end
X = temp;

clc;
for s = subs
    ns = numel(X);
    N_channels = size(X{s}{1},2);
    for Clase = 1:2
        dats = X{s}(y{s}==Clase);
        Xref = fncCutdata_(dats,0,2,fs);
        Xmi = fncCutdata_(dats,2,4,fs);
        for tr = 1:numel(dats)
            signalXref = Xref{tr};
            signalXmi = Xmi{tr};
            tic
            for ch = 1:N_channels
                sig1 = trace(cov(signalXref));
                sig2 = trace(cov(signalXmi));
                EXref(ch,:)= fncFFT(signalXref(:,ch)/sig1,fs,0);
                EXmi(ch,:) = fncFFT(signalXmi(:,ch)/sig2,fs,0);
            end
            datXmi{s}{Clase}(tr,:,:) = EXmi;
            datXref{s}{Clase}(tr,:,:) = EXref;
            fprintf(['Spec-TP: Subject ' num2str(s) ' Trial: ' num2str(tr) ' time: ' num2str(toc) '\n'])
        end
        SpecXmi{s}{Clase} = squeeze(mean(datXmi{s}{Clase},1));
    end
    SpecXref{s} = squeeze(mean([datXref{s}{1};datXref{s}{2}],1));
end
% save('BCICIV_1_TP_TR','datXmi','datXref')
% save('BCICIV1_spec22.mat','SpecXmi','SpecXref')
%% plot
load('BCICIV1_spec59.mat','SpecXmi','SpecXref')
% load('BCICIV_1_TP_TR','datXmi','datXref')
% load('electrodesBCICIV1.mat')
%
pos_ = [4,6,11,12,13,14,15,16,17,20,21,22,23,24,25,26,28,29,30,31,33,34,...
    35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,51,52,53,54,56,57,58,...
    59,60,61,62,65,66,67,68,69,70,71,76,78,85,87];
% pos_ = [4 9 10 11 12 13 15 16 17 18 19 20 21 23 24 25 26 27 31 32 33 39];
set(0,'DefaultFigureWindowStyle','docked')
L = 200; Fs = 100;
f = Fs*(0:(L/2))/L;

for s = 7
    %     for clase = 1:2
    fig=figure;
    set(fig,'name',['Subject ' num2str(s)])
    for ch = 1:59
        subplot(10,9,pos_(ch))
        %         subplot(6,7,pos_(ch))
        hold on
        %             Xrefer = [datXref{s}{1};datXref{s}{2}];
        %             for tr = 1:100
        %                 dat(tr,:) = Xrefer(tr,ch,:)/max(Xrefer(tr,ch,:));
        %             end
        %             parfor tr = 1:100
        %                 dat(tr,:) = squeeze(datXmi{s}{clase}(tr,ch,:))/max(squeeze(datXmi{s}{clase}(tr,ch,:)));
        %             end
        %             imagesc(f,1:100,squeeze(dat))
        %         plot(f,SpecXmi{s}{1}(ch,:)/max(SpecXmi{s}{1}(ch,:)),'r','linewidth',1);
        %         plot(f,SpecXmi{s}{2}(ch,:)/max(SpecXmi{s}{2}(ch,:)),'b','linewidth',1);
        %         plot(f,SpecXref{s}(ch,:)/max(SpecXref{s}(ch,:)),'g','linewidth',1);
        plot(f,SpecXmi{s}{1}(ch,:),'r','linewidth',1);
        plot(f,SpecXmi{s}{2}(ch,:),'b','linewidth',1);
        plot(f,SpecXref{s}(ch,:),'g','linewidth',1);
        hold off
        %         yticks([])
        xticks([])
        %                 title(electrodes(ch),'Interpreter','latex')
        title(num2str(ch),'Interpreter','latex')
        %         ylim([0 3*10^(-3)])
        xlim([0 50])
        %         ylim([0 100])
        drawnow
        grid on
        if ch == 22
            %             yticks([0 1])
            %                 yticks([0 1])
            xticks([0 50])
        end
    end
    
    leg1 = legend('Clase 1','Clase 2','Referencia');
    set(leg1,'Position',[0.822980277013135 0.870791312521434 0.132496511005125 0.103519665960446])
    %         suptitle(['Subject ' num2str(s) ' Class ' num2str(clase)])
    suptitle(['Subject ' num2str(s)])
    %     end
    %     subplot(10,9,88)
    %     subplot(6,7,40)
    %     colorbar
    %     axis('off')
end
