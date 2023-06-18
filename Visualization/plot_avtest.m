parfor s = 1:9
%     xpru = avtestW(s,:);
%     xtt = zeros(5,59);
%     for i=1:1
%         xtt(i,:) = xpru{:,i};
%     end
%     acc_m = mean(xtt);  % media de los folds.
%     acc_s = std(xtt);   % desviación estandar de los folds.
%     mat_w{s}(:,1) = 1:59;
%     mat_w{s}(:,2) = acc_m';        % guarda la media.
%     mat_w{s}(:,5) = acc_s';
%     mat_w{s}(:,3) = acc_m'-acc_s';
%     mat_w{s}(:,4) = acc_m'+acc_s'; % guardar la desviación estandar.
%     
    xpru = avtestTP(s,:);
    xttp = zeros(5,59);
    for i=1:1
        xttp(i,:) = xpru{:,i};
    end
    
    acc_mtp = mean(xttp);  % media de los folds.
    acc_s = std(xttp);   % desviación estandar de los folds.
    mat_tp{s}(:,1) = 1:59;
    mat_tp{s}(:,2) = acc_mtp';        % guarda la media.
    mat_tp{s}(:,5) = acc_s';
    mat_tp{s}(:,3) = acc_mtp'-acc_s';
    mat_tp{s}(:,4) = acc_mtp'+acc_s'; % guardar la desviación estandar.

    % figure
    subplot(2,2,s)
    hold on
    plot(1:59,acc_m,'r')
    plot(1:59,acc_mtp,'b')
    hold off
    grid on
    xlim([5 60])
    ylim([0.2 0.8])
    legend('wang','TP')%,'Temporal')
    if s == 1
    title(['Sujeto ' num2str(s)],'Interpreter','latex')
    elseif s == 2
        title(['Sujeto ' num2str(s)],'Interpreter','latex')
    else
       title(['Sujeto ' num2str(s+3)],'Interpreter','latex') 
    end
    xlabel('No. channels','Interpreter','latex')
    ylabel('Acc','Interpreter','latex')
end


parfor s = 1:4 max_w(s) = max(mat_w{s}(:,2)); 
    temp = find(mat_w{s}(:,2) == max_w(s));
    n_chanw(s) = temp(1);
    mm_w(s) = mat_w{s}(temp(1),2)-mat_w{s}(temp(1),3); 
end  % mejor acc.
mm_w = mm_w*100;
max_w = max_w*100;

parfor s = 1:4 matw(s,:,:) = mat_w{1,s}; end       % media del acierto por folds
matmean_w = squeeze(mean(matw,1));


parfor s = 1:4 max_tp(s) = max(mat_tp{s}(:,2)); 
    temp = find(mat_tp{s}(:,2) == max_tp(s));
    n_chantp(s) = temp(1);
    mm_tp(s) = mat_tp{s}(temp(1),2)-mat_tp{s}(temp(1),3); 
end  % mejor acc.
mm_tp = mm_tp*100;
max_tp = max_tp*100;

parfor s = 1:4 mattp(s,:,:) = mat_tp{1,s}; end       % media del acierto por folds
matmean_tp = squeeze(mean(mattp,1));