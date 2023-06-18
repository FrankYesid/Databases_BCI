% clear; clc
%% conexiones por hemisferios
clear; clc
conexiones_i = [2,3,7:9,14,15,19];   % hemisferio izquierdo.
conexiones_d = [5,6,11:13,17,18,21]; % hemisferio derecho.

c = 1;
for a = 1:22
    for b = a:22
        if a ~= b
            conexiones(c,:) = [b a];
            c = c+1;
        end
    end
end

conexion = conexiones.*0;
for i = 1:numel(conexiones_i)
    conx_i=find(conexiones(:,1)==conexiones_i(i));
    conexion(conx_i,1) = 1;
    conx_i=find(conexiones(:,2)==conexiones_i(i));
    conexion(conx_i,2) = 1;
end
ind_i = ismember(sum(conexion,2),2);

conexion = conexiones.*0;
for i = 1:numel(conexiones_d)
    conx_i=find(conexiones(:,1)==conexiones_d(i));
    conexion(conx_i,1) = 1;
    conx_i=find(conexiones(:,2)==conexiones_d(i));
    conexion(conx_i,2) = 1;
end
ind_d = ismember(sum(conexion,2),2);

conexion = conexiones.*0;
for i = 1:numel(conexiones_i)
    conx_i=find(conexiones(:,1)==conexiones_i(i));
    conexion(conx_i,1) = 1;
    conx_i=find(conexiones(:,1)==conexiones_d(i));
    conexion(conx_i,1) = 1;
    conx_i=find(conexiones(:,2)==conexiones_i(i));
    conexion(conx_i,2) = 1;
    conx_i=find(conexiones(:,2)==conexiones_d(i));
    conexion(conx_i,2) = 1;
end
ind_cr = ismember(sum(conexion,2),2);

% conexiones según el canal
temp = [0 0 0 1 0 0 0;...
                0 1 1 1 1 1 0;...
                1 1 1 1 1 1 1;...
                0 1 1 1 1 1 0;...
                0 0 1 1 1 0 0;...
                0 0 0 1 0 0 0];
plot_index = find(temp' >= 1);
n_rows = size(temp, 1);
n_cols = size(temp, 2);

counter = 1;
temp = temp';
lap = zeros(size(temp,1), size(temp,2));

% Used electrode positions instead of ones (format (1))
positions = [];
if sum(sum(temp)) ~= (sum(sum(temp>0)))
    [tmp, positions] = sort(temp(find(temp)));
    temp = temp > 0;
end

for k = 1:numel(temp)
    if temp(k) >= 1
        lap(k) = counter;
        counter = counter + 1;
    end
end

neighbors = ones(counter - 1, 4) * nan;
electrode = 0;
for (k = 1:numel(lap))
    if lap(k) ~= 0
        col = 1;
        electrode = electrode + 1;
        if (k - size(lap, 1) > 0 && lap(k - size(lap, 1)) ~= 0)  % T
            neighbors(electrode, col) = lap(k - size(lap, 1));
            col = col + 1;
        end
        if (mod(k+1, size(lap, 1)) ~= 1 && k < numel(lap) && lap(k+1) ~= 0)  % L
            neighbors(electrode, col) = lap(k+1);
            col = col + 1;
        end
        if (mod(k-1, size(lap, 1)) ~= 0 && k > 1 && lap(k-1) ~= 0)  % R
            neighbors(electrode, col) = lap(k-1);
            col = col + 1;
        end
        if (k + size(lap, 1) < numel(lap) && lap(k + size(lap, 1)) ~= 0)  % B
            neighbors(electrode, col) = lap(k + size(lap, 1));
            col = col + 1;
        end
    end
end

% conexiones respecto canales cercanos
neighbors(isnan(neighbors))= 0;
for ch = 1:size(neighbors,1)
    neig = neighbors(ch,:);    
    conexion = conexiones.*0;
    for i = 1:sum(neig>0) %numel(conexiones_i)
        conx_i=find(conexiones(:,1)==ch);
        conexion(conx_i,1) = 1;
        conx_i=find(conexiones(:,2)==ch);
        conexion(conx_i,2) = 1;
        conx_i=find(conexiones(:,1)==neig(i));
        conexion(conx_i,1) = 1;
        conx_i=find(conexiones(:,2)==neig(i));
        conexion(conx_i,2) = 1;
    end
    ind_select{ch} = ismember(sum(conexion,2),2);
end

% conexiones respecto canales 
neighbors(isnan(neighbors))= 0;
for ch = 1:size(neighbors,1)
    neig = neighbors(ch,:);
    conexion = conexiones.*0;
    for i = 1:sum(neig>0)
        for ch2 = 1:size(conexiones,1)
            if conexiones(ch2,1) == ch && conexiones(ch2,2) == neig(i) ...
                || conexiones(ch2,1) == neig(i) && conexiones(ch2,2) == ch
                conexion(ch2,1) = 1;
                conexion(ch2,2) = 1;
            end
        end
    end
    ind_select2{ch} = ismember(sum(conexion,2),2);
end


% organizar conectividades
um = 0.5;
% time
time = 0:0.04:7;
t1 = 6; t2 = 171; 

SS = 1:9;
% load('F:\BCI Competition\BCICIV_2a\Cx_funcional_BCI_ref\Cx_all.mat')
pos_ = [4,9:13,15:21,23:27,31:33,39];
mod = 3; % 1 positivos 2 negativos
coef_type = 1; % typo de grafo en clustering
for s = SS
    % name = ['Cx_wpli_BCI_2a_all_time_0_7_sub' num2str(s) '_folds_1.mat'];%['ERDs_Sub_' num2str(s) '_fold10filter_all_con_significance.mat']; name =['ERD_folds_sub_filter_graimann',num2str(s)];%
    % %         load(['D:\BCI\ERDSjunio2019\ERD_traslape\' name])
    % %         load(['D:\BCI\ERDSjunio2019\ERD_graimann\' name])
    load(['D:\Cxjulio2020\Cx_wpli_BCI_2a_all_time_0_7_sub',num2str(s),'_folds_1.mat'])
    clear Cx
    for cl = 1:numel(Cx_)
%         for mod_ = 1:2
            tem = Cx_{cl};
            tem(isnan(tem))=0;
            %             if mod_ == 1
            %                 tem(tem<0) = 0;
            %             elseif mod_ == 2
            %                 tem(tem>0) = 0;
            %                 tem = abs(tem);
            %             end
            for fr = 1:size(Cx_{cl},3)
                a = 1;
                for v = t1:t2 %1:size(Cx_{cl},4)
                    for val =1:2
                        aa  = tem(:,:,fr,v);
                        if val == 1
                            aa  = aa.*(aa>0);
                        elseif val == 2
                            aa  = -aa.*(aa<0);
                        end
                        temm= aa;%threshold_proportional(aa,um);%.*squareform(ind_select{ch});
                        temm(isnan(temm)) = 0;
                        %                     imagesc(temm); pause(0.050)
                        %                     if v == 51
                        %                         pause()
                        %                     end
                                           [C_pos,C_neg,Ctot_pos,Ctot_neg] = clustering_coef_wu_sign(temm,coef_type); % weight_conversion(temm, 'normalize')
                        %                    [C_pos,C_neg,Npos,Nneg,Kpos,Kneg] = fnc_density(temm);
                        %                    [C_pos,C_neg,vpos,vneg] = strengths_und_sign(temm);
                        %                    [C_pos,C_neg] = local_assortativity_wu_sign(temm);
                        for ch = 1:22
                            % Clustering coefficiency
%                             vall = clustering_coef_wu(weight_conversion(temm.*squareform(ind_select2{ch}), 'normalize'));
                            % Node Strength
                            vall = strengths_und(temm);
                            % Density
%                             vall = density_und(temm);
                            %eliminar NanN presentes luego de la medida
                            vall(isnan(vall)) = 0;
                            if val == 1                                
                                C_pos(ch) = vall(ch);
                            elseif val == 2
                                C_neg(ch) = vall(ch);
                            end
                        end
                    end
                    C_pos(isnan(C_pos)) = 0;
                    C_neg(isnan(C_neg)) = 0;
                    for ch = 1:22 %size(temporal,2)
                        Cxa{s}{1}{cl}{fr}(ch,a) = C_pos(ch);
                        Cxa{s}{2}{cl}{fr}(ch,a) = C_neg(ch);                        
                    end
                    a = a+1;
                end
            end
%         end
    end
end
whos
%% normalización 1: Zscore todos los canales.
for s = SS
    for mod_ = 1:2
        for cl = 1:2
            for fr = 1:17
                aaa = Cxa{s}{mod_}{cl}{fr};
                aaa(isnan(aaa))=0;
                Z1{s}{mod_}{cl}{fr} = zscore(aaa,0);
            end
        end
    end
end

%% normalización 2: min - max 
for s = SS
    for mod_ = 1:2
        maxi = max(cell2mat(cellfun(@(x) max(cell2mat(cellfun(@(x) max(x(:)),x,'UniformOutput',false))),Cxa{s}{mod_},'UniformOutput',false)));
        mini = min(cell2mat(cellfun(@(x) min(cell2mat(cellfun(@(x) min(x(:)),x,'UniformOutput',false))),Cxa{s}{mod_},'UniformOutput',false)));
        for cl = 1:2
            for fr = 1:17
                Z2{s}{mod_}{cl}{fr} = (Cxa{s}{mod_}{cl}{fr}-mini)./(maxi-mini);
            end
        end
    end
end

%% figure new
SS = 1:9;
orden = 10;
for cl = 1:2
    %         load(['F:\BCI Competition\BCICIV_2a\cx_fun_temp\WPLI_cx_new_2_all',num2str(s),num2str(cl),'_0_.mat'])
    %         subplot(6,7,pos_(ch))
    for freq = [3:5,8,10]
        figure
%         for mod_ = 1
            for ch= 1:22
                subplot(6,7,pos_(ch))
                for s = SS
                    %  subplot(6,7,pos_(ch))                    
                    hold on
                    Cx1 = Z2{s}{1}{cl}{freq}; %wpli{s,cl}.wplispctrm;
                    Cx2 = Z2{s}{2}{cl}{freq}; %wpli{s,cl}.wplispctrm;
                    %             Cx1(Cx1<0) = 0;
                    %             for freq = 1:17
                    %                 for time = 1:1750
                    %                     C_x(:,freq,time) = sum(mean(squareform(Cx(:,freq,time)),3));
                    %                 end
                    %             end
                    %             tem = squeeze(Cx(ch,:,:));
                    %             imagesc((Cx1-mini)/(maxi-mini),[0 1])
%                     rel = abs((Cx1./21)-(Cx2./(Cx1+Cx2)).*(Cx2./21));
                    rel = Cx2;
                    rel(isnan(rel))= 0;
                    %             if max(rel(:)) > maxx{s}
                    %                 maxx{s} = max(rel(:));
                    %             end
                    hold on
                    %             imagesc(0:0.0400:7,1:2:17,rel,[0,0.1])                    
                    if s == 8 || s == 3 || s == 1 || s == 9
                        c1 = [0.4660, 0.6740, 0.1880];
                        plot(time(t1:t2),medfilt1(rel(ch,:),orden),'Color',c1)
                    elseif s == 7 || s == 5 || s == 4
                        c1 = [0.9290, 0.6940, 0.1250];
                        plot(time(t1:t2),medfilt1(rel(ch,:),orden),'Color',c1)
                    elseif s == 6 || s == 2
                        c1 = [0.6350, 0.0780, 0.1840];
                        plot(time(t1:t2),medfilt1(rel(ch,:),orden),'Color',c1)
                    end
% %                     ylim([0,1])
                    %             colorbar
                    %             axis xy
                    %             axis square
                end
            end
%         end
        suptitle(['Clas ',num2str(cl),' Freq ',num2str(freq)])
%         saveas(gca,['C:\Users\frany\Desktop\figuras_ultimas',filesep,'Clas ',num2str(cl),' Freq ',num2str(freq)],'fig')
%         close
    end
end


%% figure old
set(0,'DefaultFigureWindowStyle','docked')
orden = 5; %SS = [8,3,1,9];
for cl = 1:2
    %         load(['F:\BCI Competition\BCICIV_2a\cx_fun_temp\WPLI_cx_new_2_all',num2str(s),num2str(cl),'_0_.mat'])
    %         subplot(6,7,pos_(ch))
    for freq = [3:5,8,10]
        figure
%         for mod_ = 1
            for ch= 1:22
                subplot(6,7,pos_(ch))
                for s = SS
                    %             subplot(6,7,pos_(ch))                    
                    hold on
                    Cx1 = Cxa{s}{1}{cl}{freq};%wpli{s,cl}.wplispctrm;
                    Cx2 = Cxa{s}{2}{cl}{freq};%wpli{s,cl}.wplispctrm;
                    %             Cx1(Cx1<0) = 0;
                    %             for freq = 1:17
                    %                 for time = 1:1750
                    %                     C_x(:,freq,time) = sum(mean(squareform(Cx(:,freq,time)),3));
                    %                 end
                    %             end
                    %             tem = squeeze(Cx(ch,:,:));
                    %             imagesc((Cx1-mini)/(maxi-mini),[0 1])
%                     rel = (Cx1./21)-(Cx2./(Cx1+Cx2)).*(Cx2./21);
                    rel = Cx2;
                    rel(isnan(rel))= 0;
                    %             if max(rel(:)) > maxx{s}
                    %                 maxx{s} = max(rel(:));
                    %             end
                    hold on
                    %             imagesc(0:0.0400:7,1:2:17,rel,[0,0.1])                    
                    if s == 8 || s == 3 || s == 1 || s == 9
                        c1 = [0.4660, 0.6740, 0.1880];
                        plot(time(t1:t2),medfilt1(rel(ch,:),orden),'Color',c1)
                    elseif s == 7 || s == 5 || s == 4
                        c1 = [0.9290, 0.6940, 0.1250];
                        plot(time(t1:t2),medfilt1(rel(ch,:),orden),'Color',c1)
                    elseif s == 6 || s == 2
                        c1 = [0.6350, 0.0780, 0.1840];
                        plot(time(t1:t2),medfilt1(rel(ch,:),orden),'Color',c1)
                    end
                    %             ylim([1,17])
                    %             colorbar
                    %             axis xy
                    %             axis square
                end
            end
%         end
        suptitle(['Clas ',num2str(cl),' Freq ',num2str(freq)])
%         saveas(gca,['C:\Users\frany\Desktop\figuras_ultimas',filesep,'Clas ',num2str(cl),' Freq ',num2str(freq)],'fig')
%         close
    end
end
