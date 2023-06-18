%% Figura de la sincronización por canales
% posiciones de topoplots
load('D:\Dropbox\ERD\results_ERDfc_subjects\ERDs_en topoplots\erdscolormap.mat')

SS = 1:9;
t = 0:1/250:7;
pos_ = [4,9:13,15:21,23:27,31:33,39];
for s = SS
    name = ['ERD_folds30_sub' num2str(s) '.mat'];%['ERDs_Sub_' num2str(s) '_fold10filter_all_con_significance.mat']; name =['ERD_folds_sub_filter_graimann',num2str(s)];%
    %         load(['D:\BCI\ERDSjunio2019\ERD_traslape\' name])
    %         load(['D:\BCI\ERDSjunio2019\ERD_graimann\' name])
    load(['I:\ERD\' name])
    for cl = 1:2
        for ch= 1:22
            subplot(6,7,pos_(ch))
            imagesc(ERDsfilt_{1}{cl}{ch},[-1, 1.5])
            axis xy
            colormap(erdcolormap)
            colorbar
        end
        saveas(gca,['D:\ERDs_sub',num2str(s),'_clas_',num2str(cl)],'png')
        %         close
    end
    clear ERDsfilt_
end

%% Figura de la conectividad por canales strength aplicando una normalización
SS = [8,2]; um = 0.3;
% load('F:\BCI Competition\BCICIV_2a\Cx_funcional_BCI_ref\Cx_all.mat')
pos_ = [4,9:13,15:21,23:27,31:33,39];
mod = 1; % 1 positivos 2 negativos
for s = SS
    % name = ['Cx_wpli_BCI_2a_all_time_0_7_sub' num2str(s) '_folds_1.mat'];%['ERDs_Sub_' num2str(s) '_fold10filter_all_con_significance.mat']; name =['ERD_folds_sub_filter_graimann',num2str(s)];%
    % %         load(['D:\BCI\ERDSjunio2019\ERD_traslape\' name])
    % %         load(['D:\BCI\ERDSjunio2019\ERD_graimann\' name])
    load(['D:\Cxjulio2020\Cx_wpli_BCI_2a_all_time_0_7_sub',num2str(s),'_folds_1.mat'])
    clear Cx
    for cl = 1:numel(Cx_)
        tem = Cx_{cl};
        tem(isnan(tem))=0;
        %         for mod_ = 1:2
        if mod == 1
            tem(tem<0) = 0;
        elseif mod == 2
            tem(tem>0) = 0;
            tem = abs(tem);
        end
        for fr = 1:size(Cx_{cl},3)
            for v = 1:size(Cx_{cl},4)
                val = threshold_proportional(tem(:,:,fr,v),um);
                temporal = strengths_und(val);
                for ch = 1:size(temporal,2)
                    Cx{cl}{ch}(fr,v) = temporal(ch);
                end
            end
        end
        %         end
    end
    %     for mod_ = 1:2
    maxi = max(cell2mat(cellfun(@(x) max(cell2mat(cellfun(@(y) max(y(:)),x,'UniformOutput',false))),Cx,'UniformOutput',false)));
    mini = 0;
    %     end
    %     maxx{s} = 0;
    %     maxx{2} = 0.0841;
    %     maxx{8} = 0.0992;
    for ch= 1:22%cl = 1:2
        %         figure
        %         load(['F:\BCI Competition\BCICIV_2a\cx_fun_temp\WPLI_cx_new_2_all',num2str(s),num2str(cl),'_0_.mat'])
        subplot(6,7,pos_(ch))
        for cl = 1:2%ch= 1:22
            %             subplot(6,7,pos_(ch))
            hold on
            Cx1 = Cx{cl}{ch};%wpli{s,cl}.wplispctrm;
            Cx2 = Cx{cl}{ch};%wpli{s,cl}.wplispctrm;
            %             Cx1(Cx1<0) = 0;
            %             for freq = 1:17
            %                 for time = 1:1750
            %                     C_x(:,freq,time) = sum(mean(squareform(Cx(:,freq,time)),3));
            %                 end
            %             end
            %             tem = squeeze(Cx(ch,:,:));
            %             imagesc((Cx1-mini)/(maxi-mini),[0 1])
            rel = Cx1./maxi;
            rel(isnan(rel))= 0;
            %             if max(rel(:)) > maxx{s}
            %                 maxx{s} = max(rel(:));
            %             end
            hold on
            %             imagesc(0:0.0400:7,1:2:17,rel,[0,1])
            plot(0:0.0400:7,rel(3,:))
            %             ylim([1,17])
            %             colorbar
            %             axis xy
            %             axis square
        end
        %         legend('Class 1','Class 2')
        %         if mod == 1
        %             suptitle(['Sujeto ',num2str(s),' cl ',num2str(cl),' pos'])
        %         elseif mod == 2
        %             suptitle(['Sujeto ',num2str(s),' cl ',num2str(cl),' neg'])
        %         end
        %         saveas(gca,['D:\Dropbox\[3.1] MI ERDs Mask\figures\BCI_CIV_2a\comparación_ERD_Cx\wpli_ult_sub',num2str(s),'_clas_',num2str(cl),'_',num2str(mod)],'png')
        %         saveas(gca,['D:\Dropbox\[3.1] MI ERDs Mask\figures\BCI_CIV_2a\comparación_ERD_Cx\wpli_ult_sub',num2str(s),'_clas_',num2str(cl),'_',num2str(mod)],'fig')
        %         close
    end
    suptitle(['Sujeto ',num2str(s),' pos'])
    saveas(gca,['D:\Dropbox\[3.1] MI ERDs Mask\figures\BCI_CIV_2a\comparación_ERD_Cx\wpli_ult_plot_sub',num2str(s),'_clas_',num2str(cl),'_',num2str(mod)],'png')
end


%% Figura de la conectividad por canales strength aplicando una normalización según robinov y sporns
SS = [8,2];
% load('F:\BCI Competition\BCICIV_2a\Cx_funcional_BCI_ref\Cx_all.mat')
pos_ = [4,9:13,15:21,23:27,31:33,39];
mod = 3; % 1 positivos 2 negativos
for s = SS
    % name = ['Cx_wpli_BCI_2a_all_time_0_7_sub' num2str(s) '_folds_1.mat'];%['ERDs_Sub_' num2str(s) '_fold10filter_all_con_significance.mat']; name =['ERD_folds_sub_filter_graimann',num2str(s)];%
    % %         load(['D:\BCI\ERDSjunio2019\ERD_traslape\' name])
    % %         load(['D:\BCI\ERDSjunio2019\ERD_graimann\' name])
    load(['D:\Cxjulio2020\Cx_wpli_BCI_2a_all_time_0_7_sub',num2str(s),'_folds_1.mat'])
    clear Cx
    for cl = 1:numel(Cx_)
        tem = Cx_{cl};
        tem(isnan(tem))=0;
        for mod_ = 1:2
            if mod_ == 1
                tem(tem<0) = 0;
            elseif mod_ == 2
                tem(tem>0) = 0;
                tem = abs(tem);
            end
            for fr = 1:size(Cx_{cl},3)
                for v = 1:size(Cx_{cl},4)
                    temporal = strengths_und(tem(:,:,fr,v));
                    for ch = 1:size(temporal,2)
                        Cx{mod_}{cl}{ch}(fr,v) = temporal(ch);
                    end
                end
            end
        end
    end
    for mod_ = 1:2
        maxi = max(cell2mat(cellfun(@(x) max(cell2mat(cellfun(@(y) max(y(:)),x,'UniformOutput',false))),Cx{mod_},'UniformOutput',false)));
        mini = 0;
    end
    %     maxx{s} = 0;
    maxx{2} = 0.0841;
    maxx{8} = 0.0992;
    for ch= 1:22
        %         figure
        %         load(['F:\BCI Competition\BCICIV_2a\cx_fun_temp\WPLI_cx_new_2_all',num2str(s),num2str(cl),'_0_.mat'])
        subplot(6,7,pos_(ch))
        for cl = 1:2
            hold on
            Cx1 = Cx{1}{cl}{ch};%wpli{s,cl}.wplispctrm;
            Cx2 = Cx{2}{cl}{ch};%wpli{s,cl}.wplispctrm;
            %             Cx1(Cx1<0) = 0;
            %             for freq = 1:17
            %                 for time = 1:1750
            %                     C_x(:,freq,time) = sum(mean(squareform(Cx(:,freq,time)),3));
            %                 end
            %             end
            %             tem = squeeze(Cx(ch,:,:));
            %             imagesc((Cx1-mini)/(maxi-mini),[0 1])
            rel = (Cx1./63)-(Cx2./(Cx1+Cx2)).*(Cx2./63);
            rel(isnan(rel))= 0;
            %             if max(rel(:)) > maxx{s}
            %                 maxx{s} = max(rel(:));
            %             end
            %             imagesc(0:0.0400:7,1:2:17,rel./maxx{s},[0,1])
            plot(0:0.0400:7,rel(3,:)./maxx{s})
            %             ylim([1,17])
            %             colorbar
            %             axis xy
            %             axis square
        end
        %         if mod == 1
        %             suptitle(['Sujeto ',num2str(s),' cl ',num2str(cl),' pos'])
        %         elseif mod == 2
        %             suptitle(['Sujeto ',num2str(s),' cl ',num2str(cl),' neg'])
        %         elseif mod == 3
        %             suptitle(['Sujeto ',num2str(s),' cl ',num2str(cl),' norm'])
        %         end
        %         saveas(gca,['D:\Dropbox\[3.1] MI ERDs Mask\figures\BCI_CIV_2a\comparación_ERD_Cx\wpli_ult2_sub',num2str(s),'_clas_',num2str(cl),'_',num2str(mod)],'png')
        %         saveas(gca,['D:\Dropbox\[3.1] MI ERDs Mask\figures\BCI_CIV_2a\comparación_ERD_Cx\wpli_ult2_sub',num2str(s),'_clas_',num2str(cl),'_',num2str(mod)],'fig')
        %         close
    end    
    suptitle(['Sujeto ',num2str(s),' pos'])
    saveas(gca,['D:\Dropbox\[3.1] MI ERDs Mask\figures\BCI_CIV_2a\comparación_ERD_Cx\wpli_ult2_plot_sub',num2str(s),'_clas_',num2str(cl),'_',num2str(mod)],'png')
end
