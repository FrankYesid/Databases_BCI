function MyTopo_fun(rel,sel,pos,label,lims,cur,ticks,t_e,database,HeadModel,type,tex,Color)
% %% Function to graph the different topoplots of the database. %%
% MyTopo_fun(Y,sel,pos,label,lims,cur,ticks,t_e,database,HeadModel)
% rel       -   Vector of weights 1xN.
% pos       -   Positions vector Nx2.
% label     -   tag cell Nx1.
% lims      -   vector colorbar limits [min max].
% cur       -   boolean enables contour lines 1,0.
% ticks     -   ability to plot Boolean labels 1,0.
% t_e       -   size of the electrodes.
% database  -   selected database.
% HeadModel -   model of the human head.
% type      -   type of topoplot.
% tex       -   ability to plot labels of channels 1,0.
% Example:
%       database = 2;  Database of BCICIV_1
%       rel = [1:22]; rel(1:22)=22;
%       load('BCICIV_1\electrodesBCICIV1.mat') % se encuentra la ubicacion de los electrodos y sus nombres.
%       load('HeadModel.mat')   % model of the head.
%       M1.xy = pos(:,1:2);     % posicion de los canales.
%       M1.lab = electrodes;    % nombre de los canales.
%       t_e = 200;
%       type = 1;
%       tex = 0;
%       MyTopo_fun(rel,sel,M1.xy,M1.lab,[min(rel) max(rel)],0,0,t_e,database,HeadModel,type,tex)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Topoplot: F.Y. Zapata-castaño
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch database
    
    case 1
        switch type
            
            case 1
                for i=1:2
                    pos(:,i) = 1.*((pos(:,i)-min(pos(:,i)))/range(pos(:,i))-0.5);
                end
                xc = HeadModel(1,:);
                yc = HeadModel(2,:);
                tam = 1;
                %                 plot(tam*xc,tam*yc,'k','LineWidth',4)
                hold on
                % Topoplot
                x = pos(:,1);
                y = pos(:,2);
                plot(0.99.*xc,0.99.*yc,'Color',[130,130,130]/255,'LineWidth',6)
                
                scatter(0.99.*x(sel),0.99.*y(sel),t_e,'b','filled')
                if tex == 1
                    text(pos(sel,1)-0.02,pos(sel,2)-0.04,label(sel),'Interpreter',...
                        'latex','ButtonDownFcn',{@lineCallback,database},'Color','black',...
                        'FontSize',10);
                end
                hold off
                
            case 2
                for i=1:2
                    pos(:,i) = 1.*((pos(:,i)-min(pos(:,i)))/range(pos(:,i))-0.5);
                end
                xc = HeadModel(1,:);
                yc = HeadModel(2,:);
                tam = 1;
                plot(0.99.*xc,0.99.*yc,'Color',[130,130,130]/255,'LineWidth',6)
                hold on
                % Topoplot
                x = pos(:,1);
                y = pos(:,2);
                scatter(0.99.*x(sel),0.99.*y(sel),t_e*rel,'b','filled')
                if tex == 1
                    text(pos(sel,1)-0.02,pos(sel,2)-0.04,label(sel),'Interpreter',...
                        'latex','ButtonDownFcn',{@lineCallback,database},'Color','black',...
                        'FontSize',10);
                end
                hold off
                
            case 3
                %                 for i=1:2
                %                     pos(:,i) = 1.*((pos(:,i)-min(pos(:,i)))/range(pos(:,i))-0.5);
                %                 end
                %                 xc = HeadModel(1,:);
                %                 yc = HeadModel(2,:);
                %                 tam = 1;
                %                 plot(tam*xc,tam*yc,'k','LineWidth',4)
                %                 hold on
                %                 % Topoplot
                %                 x = pos(:,1);
                %                 y = pos(:,2);
                fprintf('Falta \n')
            case 4
                
                for i=1:2
                    pos(:,i) = 1.*((pos(:,i)-min(pos(:,i)))/range(pos(:,i))-0.5);
                end
                xc = HeadModel(1,:);
                yc = HeadModel(2,:);
                tam = 1;
                % plot(tam*xc,tam*yc,'k','LineWidth',4)
                hold on
                % Topoplot
                x = pos(:,1);
                y = pos(:,2);
                GS = 1000;
                xi = linspace(min(x)-0.05, max(x)+0.05, GS);       % x-axis for interpolation (row vector)
                yi = linspace(min(y)-0.05, max(y)+0.05, GS);       % y-axis for interpolation (row vector)
                [Xi, Yi, Zi] = griddata([x' xc], [y' yc], [rel(:);zeros(numel(xc),1)], xi', yi,'v4'); % interpolate the topographic data
                %% Creating data mask
                [TH R] = cart2pol(Xi,Yi);
                Zi(R>0.5) = NaN;
                deltax = xi(2)-xi(1); % length of grid entry
                deltay = yi(2)-yi(1); % length of grid entry
                h = surf(Xi-deltax/2, Yi-deltay/2+0.004, zeros(size(Zi)), Zi,'EdgeColor', 'none', 'FaceColor', 'flat');hold on
                shading interp
                plot(0.99.*xc,0.99.*yc,'Color',[130,130,130]/255,'LineWidth',6)
                
                scatter(0.99.*x(sel),0.99.*y(sel),t_e,'b','filled')
                if tex == 1
                    text(pos(sel,1)-0.02,pos(sel,2)-0.04,label(sel),'Interpreter',...
                        'latex','ButtonDownFcn',{@lineCallback,database},'Color','black',...
                        'FontSize',10);
                end
                hold off
        end
        
    case 2
        switch type
            
            case 1
                for i=1:2
                    pos(:,i) = 0.8.*((pos(:,i)-min(pos(:,i)))/range(pos(:,i))-0.5);
                end
                xc = HeadModel(1,:);
                yc = HeadModel(2,:);
                tam = 1;
                % plot(tam*xc,tam*yc,'k','LineWidth',4)
                hold on
                % Topoplot
                x = pos(:,1);
                y = pos(:,2);
                tmp = [x,y,x*0]*rotz(2);
                x = tmp(:,1);
                y = tmp(:,2);
                pos(:,1) = x;
                pos(:,2) = y;
                scatter(x(sel),y(sel)-0.08,t_e,'b','filled');
                plot(tam*xc,tam*yc,'Color',[130,130,130]/255,'LineWidth',6)
                if tex == 1
                    text(pos(sel,1)-0.02,pos(sel,2)-0.12,label(sel),'Interpreter',...
                        'latex','ButtonDownFcn',{@lineCallback,database},'Color','black',...
                        'FontSize',10);
                end
                hold off
                
            case 2
                for i=1:2
                    pos(:,i) = 0.8.*((pos(:,i)-min(pos(:,i)))/range(pos(:,i))-0.5);
                end
                xc = HeadModel(1,:);
                yc = HeadModel(2,:);
                tam = 1;
                % plot(tam*xc,tam*yc,'k','LineWidth',4)
                hold on
                % Topoplot
                x = pos(:,1);
                y = pos(:,2);
                tmp = [x,y,x*0]*rotz(2);
                x = tmp(:,1);
                y = tmp(:,2);
                pos(:,1) = x;
                pos(:,2) = y;
                scatter(x(sel),y(sel)-0.08,t_e*rel,'b','filled');
                plot(tam*xc,tam*yc,'Color',[130,130,130]/255,'LineWidth',6)
                if tex == 1
                    text(pos(sel,1)-0.02,pos(sel,2)-0.12,label(sel),'Interpreter',...
                        'latex','ButtonDownFcn',{@lineCallback,database},'Color','black',...
                        'FontSize',10);
                end
                hold off
                
            case 3
                %                 for i=1:2
                %                     pos(:,i) = 1.*((pos(:,i)-min(pos(:,i)))/range(pos(:,i))-0.5);
                %                 end
                %                 xc = HeadModel(1,:);
                %                 yc = HeadModel(2,:);
                %                 tam = 1;
                %                 plot(tam*xc,tam*yc,'k','LineWidth',4)
                %                 hold on
                %                 % Topoplot
                %                 x = pos(:,1);
                %                 y = pos(:,2);
                fprintf('Falta \n')
            case 4
                
                for i=1:2
                    pos(:,i) = 0.9.*((pos(:,i)-min(pos(:,i)))/range(pos(:,i))-0.5);
                end
                xc = HeadModel(1,:);
                yc = HeadModel(2,:);
                tam = 1;
                %                 plot(tam*xc,tam*yc,'k','LineWidth',4)
                hold on
                % Topoplot
                x = pos(:,1);
                y = pos(:,2);
                tmp = [x,y,x*0]*rotz(2);
                x = tmp(:,1);
                y = tmp(:,2);
                pos(:,1) = x;
                pos(:,2) = y;
                GS = 1000;
                xi = linspace(min(x)-0.05, max(x)+0.05, GS);       % x-axis for interpolation (row vector)
                yi = linspace(min(y)-0.05, max(y)+0.05, GS);       % y-axis for interpolation (row vector)
                [Xi, Yi, Zi] = griddata([x' xc], [y' yc], [rel(:);zeros(numel(xc),1)], xi', yi,'v4'); % interpolate the topographic data
                %% Creating data mask
                [TH R] = cart2pol(Xi,Yi);
                Zi(R>0.5) = NaN;
                deltax = xi(2)-xi(1); % length of grid entry
                deltay = yi(2)-yi(1); % length of grid entry
                 %% plot acc y kappa
                figure; x1 = 0.02;     y1 = 0.02;     w = 0.95;      h = 0.95;
                subplot('position',[x1,y1,w,h]);
                surf(Xi-deltax/2, Yi-deltay/2+0.004, zeros(size(Zi)), Zi,'EdgeColor', 'none', 'FaceColor', 'flat');hold on
                
                caxis([lims(1) lims(2)])
                colormap(Color)
                shading interp
                plot(0.999.*xc,0.999.*yc,'Color',[130,130,130]/255,'LineWidth',2) % cambié el tamaño del contorno
                view([0,90])
                axis('square')
                scatter(0.999.*x(sel),0.999.*y(sel),t_e,'b','filled')
                if tex == 1
                    text(pos(sel,1)-0.02,pos(sel,2)-0.04,label(sel),'Interpreter',...
                        'latex','ButtonDownFcn',{@lineCallback,database},'Color','black',...
                        'FontSize',10);
                end
                hold off
                axis off
        end
        
    case 3
        switch type
            
            case 1
                for i=1:2
                    pos(:,i) = 1.*((pos(:,i)-min(pos(:,i)))/range(pos(:,i))-0.5);
                end
                xc = HeadModel(1,:);
                yc = HeadModel(2,:);
                tam = 1;
                %                 plot(tam*xc,tam*yc,'k','LineWidth',4)
                hold on
                % Topoplot
                x = pos(:,1);
                y = pos(:,2);
                plot(0.99.*xc,0.99.*yc,'Color',[130,130,130]/255,'LineWidth',6)
                
                scatter(x(sel),y(sel),t_e,'b','filled')
                if tex == 1
                    text(pos(sel,1)-0.02,pos(sel,2)-0.04,label(sel),'Interpreter',...
                        'latex','ButtonDownFcn',{@lineCallback,database},'Color','black',...
                        'FontSize',10);
                end
                hold off
                
            case 2
                for i=1:2
                    pos(:,i) = 1.*((pos(:,i)-min(pos(:,i)))/range(pos(:,i))-0.5);
                end
                xc = HeadModel(1,:);
                yc = HeadModel(2,:);
                tam = 1;
                plot(0.99.*xc,0.99.*yc,'Color',[130,130,130]/255,'LineWidth',6)
                hold on
                % Topoplot
                x = pos(:,1);
                y = pos(:,2);
                scatter(0.99.*x(sel),0.99.*y(sel),t_e,'b','filled')
                if tex == 1
                    text(pos(sel,1)-0.02,pos(sel,2)-0.04,label(sel),'Interpreter',...
                        'latex','ButtonDownFcn',{@lineCallback,database},'Color','black',...
                        'FontSize',10);
                end
                hold off
                
            case 3
                %                 for i=1:2
                %                     pos(:,i) = 1.*((pos(:,i)-min(pos(:,i)))/range(pos(:,i))-0.5);
                %                 end
                %                 xc = HeadModel(1,:);
                %                 yc = HeadModel(2,:);
                %                 tam = 1;
                %                 plot(tam*xc,tam*yc,'k','LineWidth',4)
                %                 hold on
                %                 % Topoplot
                %                 x = pos(:,1);
                %                 y = pos(:,2);
                fprintf('Falta \n')
            case 4
                
                for i=1:2
                    pos(:,i) = 1.*((pos(:,i)-min(pos(:,i)))/range(pos(:,i))-0.5);
                end
                xc = HeadModel(1,:);
                yc = HeadModel(2,:);
                tam = 1;
                %                 plot(tam*xc,tam*yc,'k','LineWidth',4)
                hold on
                % Topoplot
                x = pos(:,1);
                y = pos(:,2);
                GS = 1000;
                xi = linspace(min(x)-0.05, max(x)+0.05, GS);       % x-axis for interpolation (row vector)
                yi = linspace(min(y)-0.05, max(y)+0.05, GS);       % y-axis for interpolation (row vector)
                [Xi, Yi, Zi] = griddata([x' xc], [y' yc], [rel(:);zeros(numel(xc),1)], xi', yi,'v4'); % interpolate the topographic data
                %% Creating data mask
                [TH R] = cart2pol(Xi,Yi);
                Zi(R>0.5) = NaN;
                deltax = xi(2)-xi(1); % length of grid entry
                deltay = yi(2)-yi(1); % length of grid entry
                h = surf(Xi-deltax/2, Yi-deltay/2+0.004, zeros(size(Zi)), Zi,'EdgeColor', 'none', 'FaceColor', 'flat');hold on
                shading interp
                plot(0.99.*xc,0.99.*yc,'Color',[130,130,130]/255,'LineWidth',6)
                
                scatter(0.99.*x(sel),0.99.*y(sel),t_e,'b','filled')
                if tex == 1
                    text(pos(sel,1)-0.02,pos(sel,2)-0.04,label(sel),'Interpreter',...
                        'latex','ButtonDownFcn',{@lineCallback,database},'Color','black',...
                        'FontSize',10);
                end
                hold off
        end
        
    case 4
        switch type
            
            case 1
                for i=1:2
                    pos(:,i) = 1.*((pos(:,i)-min(pos(:,i)))/range(pos(:,i))-0.5);
                end
                xc = HeadModel(1,:);
                yc = HeadModel(2,:);
                tam = 1;
                %                 plot(tam*xc,tam*yc,'k','LineWidth',4)
                hold on
                % Topoplot
                x = pos(:,1);
                y = pos(:,2);
                plot(0.99.*xc,0.99.*yc,'Color',[130,130,130]/255,'LineWidth',6)
                
                scatter(x(sel),y(sel),t_e,'b','filled')
                if tex == 1
                    text(pos(sel,1)-0.02,pos(sel,2)-0.04,label(sel),'Interpreter',...
                        'latex','ButtonDownFcn',{@lineCallback,database},'Color','black',...
                        'FontSize',10);
                end
                hold off
                
            case 2
                for i=1:2
                    pos(:,i) = 1.*((pos(:,i)-min(pos(:,i)))/range(pos(:,i))-0.5);
                end
                xc = HeadModel(1,:);
                yc = HeadModel(2,:);
                tam = 1;
                plot(0.99.*xc,0.99.*yc,'Color',[130,130,130]/255,'LineWidth',6)
                hold on
                % Topoplot
                x = pos(:,1);
                y = pos(:,2);
                scatter(0.99.*x(sel),0.99.*y(sel),t_e,'b','filled')
                if tex == 1
                    text(pos(sel,1)-0.02,pos(sel,2)-0.04,label(sel),'Interpreter',...
                        'latex','ButtonDownFcn',{@lineCallback,database},'Color','black',...
                        'FontSize',10);
                end
                hold off
                
            case 3
                %                 for i=1:2
                %                     pos(:,i) = 1.*((pos(:,i)-min(pos(:,i)))/range(pos(:,i))-0.5);
                %                 end
                %                 xc = HeadModel(1,:);
                %                 yc = HeadModel(2,:);
                %                 tam = 1;
                %                 plot(tam*xc,tam*yc,'k','LineWidth',4)
                %                 hold on
                %                 % Topoplot
                %                 x = pos(:,1);
                %                 y = pos(:,2);
                fprintf('Falta \n')
            case 4
                
                for i=1:2
                    pos(:,i) = 1.*((pos(:,i)-min(pos(:,i)))/range(pos(:,i))-0.5);
                end
                xc = HeadModel(1,:);
                yc = HeadModel(2,:);
                tam = 1;
                %                 plot(tam*xc,tam*yc,'k','LineWidth',4)
                hold on
                % Topoplot
                x = pos(:,1);
                y = pos(:,2);
                GS = 1000;
                xi = linspace(min(x)-0.05, max(x)+0.05, GS);       % x-axis for interpolation (row vector)
                yi = linspace(min(y)-0.05, max(y)+0.05, GS);       % y-axis for interpolation (row vector)
                [Xi, Yi, Zi] = griddata([x' xc], [y' yc], [rel(:);zeros(numel(xc),1)], xi', yi,'v4'); % interpolate the topographic data
                %% Creating data mask
                [TH R] = cart2pol(Xi,Yi);
                Zi(R>0.5) = NaN;
                deltax = xi(2)-xi(1); % length of grid entry
                deltay = yi(2)-yi(1); % length of grid entry
                h = surf(Xi-deltax/2, Yi-deltay/2+0.004, zeros(size(Zi)), Zi,'EdgeColor', 'none', 'FaceColor', 'flat');hold on
                shading interp
                plot(0.99.*xc,0.99.*yc,'Color',[130,130,130]/255,'LineWidth',6)
                
                scatter(0.99.*x(sel),0.99.*y(sel),t_e,'b','filled')
                if tex == 1
                    text(pos(sel,1)-0.02,pos(sel,2)-0.04,label(sel),'Interpreter',...
                        'latex','ButtonDownFcn',{@lineCallback,database},'Color','black',...
                        'FontSize',10);
                end
                hold off
        end
end




