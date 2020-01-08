function plotSyncValuesVsBic_paper_only10uM(method, methodNames)



%% 2) Plot
hf=figure(1);
hf.Name='Figure: Results';
hf.Units = 'centimeter';
H = 18;
W = 12; %size(method,2)*10/2;
hf.Position = [0,0,W,H]; % change size of displayed window
hf.PaperPosition=[0,0,W,H]; % change size of file when "print" is used


ROW=2;                  % sync-curve, ttest
COL=1;     % all methods into one plot
x = 1:size(method,2);
colorRef = [0.8 0.8 0.8];
colorBic = [0 0 0];
color1 = [0 0.5 1]; % color for 5% significance line
color2 = [0 1 0.5]; % color for 10% significance line

for m=1:size(method,2)
    
    
    %% 1) Sync-Curve
    row=1;
    hs(row)=subplot(ROW,COL,row);
    REF = method(m).REF;
    BIC = method(m).BIC;
    x = [1+(m-1)*3, 2+(m-1)*3];
    x_middle(m) = x(1)+ (x(2)-x(1))/2;
    for i=1:size(REF,2)
        plot(x,[REF(i), BIC(i)],'-','Color',colorRef); hold on
        hp(1)=plot(x(1),REF(i),'o','Color',colorBic); hold on
        hp(2)=plot(x(2),BIC(i),'o','Color',colorBic); hold on
        
        hp(2).MarkerFaceColor=[0.5 0.5 0.5];
    end
    hs(row).YLim=[0 1];
    for i=1:length(hs(row).YTick)
        hs(row).YTickLabel{i} = num2str(hs(row).YTick(i),'%.1f');
    end
    if m==1; hs(row).YLabel.String={'Synchrony'}; end
    method(m).name = strrep(method(m).name,'_','-');
    method(m).name = strrep(method(m).name,'Sync-','');
    if m==size(method,2)
        hs(row).XTick = x_middle;
        hs(row).XTickLabel = methodNames;
        hs(row).XTickLabelRotation=45;
        legend(hp,'0 µM BIC','10 µM BIC','LOCATION','SouthEast')
    end
    hs(row).YLim=[0 1];
    
    
    
    
    %% 3) ttest
    row=row+1;
    hs(row)=subplot(ROW,COL,row);
    tmp=method(m).pp;
    hp=plot(x_middle(m),tmp,'kx'); hold on
    hp.MarkerSize=8;
    hp.LineWidth=1.0;
    hs(row).YScale='log';
    hs(row).YLim=[1e-4 1];
    hs(row).YTick=[1e-4 1e-2 1];
    if m==1; hs(row).YLabel.String='p (t-test)'; end
    
    if m==size(method,2)
        hp=plot([hs(row).XLim(1),hs(row).XLim(2)],[0.05 0.05],'--','Color',color1);
        hp=plot([hs(row).XLim(1),hs(row).XLim(2)],[0.01 0.01],'--','Color',color2);
        
        text(hs(row).XLim(2),0.05,'5% (*)','FontSize',8,'Color',color1)
        text(hs(row).XLim(2),0.01,'1% (**)','FontSize',8,'Color',color2)
        hs(row).XTick = x_middle;
        hs(row).XTickLabel = methodNames;
        hs(row).XTickLabelRotation=45;
        
        hs(row).YTick=[10^-4,10^-3,10^-2,10^-1,10^0];
        hs(row).YGrid = 'on';
        
    end
    
    
end

end