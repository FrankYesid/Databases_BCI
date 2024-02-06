clear
clc

% load('C:\Users\frany\Dropbox\Event-related\Codes\TP\Matlab_con_csp\BCICIV_1.mat')
load('G:\Dropbox\ERD\Codes\TP\Matlab_con_csp\BCICIV_1.mat')
fs = 100;
seg_start = 50;
seg_end = 250;
subs = [2,7];

% load('BCICIV_2a.mat')
% seg_start = 250;
% seg_end = 750;
% subs = [2,3];

Freq = [8,30];
temp=cell(numel(subs),1);
for s = subs
    temp{s} = fcnfiltband(X{s},fs,Freq,8);
end
X = temp;

%% modo 1
for s = subs
    ref = cellfun(@(x)(x(1:seg_start,:)'),X{s},'UniformOutput',false);
    signalcl1 = X{s}(y{s}==1);
    sigref1 = cellfun(@(x)(x(1:seg_start,:)'),signalcl1,'UniformOutput',false);
    signalcl1 = cellfun(@(x)(x(seg_start+1:seg_end,:)'),signalcl1,'UniformOutput',false);
    signalcl2 = X{s}(y{s}==2);
    sigref2 = cellfun(@(x)(x(1:seg_start,:)'),signalcl2,'UniformOutput',false);
    signalcl2 = cellfun(@(x)(x(seg_start+1:seg_end,:)'),signalcl2,'UniformOutput',false);
    W{s} = fncCSPcohen( ref,[signalcl1;signalcl2], 2, 3); % Finding
    components{s} = cellfun(@(x)(x*W{s}'/trace(cov(x))),X{s},'UniformOutput',false);
end
% save('Components_22.mat','components')
%% modo 2
for s = subs
    signalcl1 = X{s}(y{s}==1);
    signalcl2 = X{s}(y{s}==2);
    if numel(signalcl1) > numel(signalcl2)
        te = signalcl1;
        clear signalcl1
        signalcl1 = te(1:numel(signalcl2));
    elseif numel(signalcl1) < numel(signalcl2)
        te = signalcl2;
        clear signalcl2
        signalcl2 = te(1:numel(signalcl1));
    end
    % clase 1
    signalcl1 = cellfun(@(x)(x(200+1:seg_end,:)),signalcl1,'UniformOutput',false);
    sig1 = ord(signalcl1);
    sigref1 = cellfun(@(x)(x(1:seg_start,:)'),signalcl1,'UniformOutput',false);
    sig1r = ord(sigref1);
    % clase 2
    signalcl2 = cellfun(@(x)(x(200+1:seg_end,:)),signalcl2,'UniformOutput',false);
    sig2 = ord(signalcl2);
    sigref2 = cellfun(@(x)(x(1:seg_start,:)'),signalcl2,'UniformOutput',false);
    sig2r = ord(sigref2);
    for ch = 1:numel(sig1)
        W{s,ch} = fncCSPcohen(sig1(ch),sig2(ch),2,3);
%         components{s,ch} = cellfun(@(x)(x*W{s,ch}'/trace(cov(x))),X{s},'UniformOutput',false);
    end
end
%%
for s = subs
    ns = numel(X);
    %     N_channels = size(X{s}{1},2);
    for Clase = 1:2
        dats = components{s}(y{s}==Clase);
        Xref =  cellfun(@(x)(x(1:seg_start,:)),dats,'UniformOutput',false);
        Xmi =  cellfun(@(x)(x(seg_start+1:seg_end,:)),dats,'UniformOutput',false);
        for tr = 1:numel(dats)
            signalXref = Xref{tr};
            signalXmi = Xmi{tr};
            tic
            for ch = 1%:N_channels
                EXref(ch,:)= fncFFT(signalXref,fs,0);
                EXmi(ch,:) = fncFFT(signalXmi,fs,0);
            end
            datXmi{s}{Clase}(tr,:,:) = EXmi;
            datXref{s}{Clase}(tr,:,:) = EXref;
            fprintf(['Spec-TP: Subject ' num2str(s) ' Trial: ' num2str(tr) ' time: ' num2str(toc) '\n'])
        end
        SpecXmi{s}{Clase} = squeeze(mean(datXmi{s}{Clase},1));
    end
end

for s = subs
    SpecXref{s} = squeeze(mean([datXref{s}{1};datXref{s}{2}],1));
end

%%
set(0,'DefaultFigureWindowStyle','docked')
L = seg_start; %Fs = 250;
f = fs*(0:(L/2))/L;

for s = subs
    fig=figure;
    set(fig,'name',['Subject ' num2str(s)])
    hold on
    plot(f,SpecXmi{s}{1},'r','linewidth',1);
    plot(f,SpecXmi{s}{2},'b','linewidth',1);
    plot(f,SpecXref{s},'g','linewidth',1);
    hold off
    xlim([0 50])
    drawnow
    grid on
    leg1 = legend('Clase 1','Clase 2','Referencia');
    set(leg1,'Position',[0.822980277013135 0.870791312521434 0.132496511005125 0.103519665960446])
    suptitle(['Subject ' num2str(s)])
end

