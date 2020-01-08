function plotSyncValuesVsBic_paper_only10uM_VarTh(th, methodNames)

for t=1:size(th,2)
    
    %% 2) Plot
    hf=figure(1);
    hf.Name='Figure: Results';
    hf.Units = 'centimeter';
    H = 12;
    W = 12; %size(method,2)*10/2;
    hf.Position = [0,0,W,H]; % change size of displayed window
    hf.PaperPosition=[0,0,W,H]; % change size of file when "print" is used
    
    
    ROW=1;                  % sync-curve, ttest, wsr-test
    COL=1;%size(method,2);     % all methods into one plot
    x = 1:size(th(t).method,2);
    colorRef = [0.8 0.8 0.8];
    colorBic = [0 0 0];
    color1 = [0 0.5 1]; % color for 5% significance line
    color2 = [0 1 0.5]; % color for 10% significance line
    
    for m=1:size(th(t).method,2)
        
        [color,line,marker]=getPlotStyle(m);
        
        % x
        x = 1:1:size(th,2);
        
        % y
        for t2=1:size(th,2)
            pp(t2)=th(t2).method(m).pp;
            ppf(t2)=th(t2).method(m).ppf;
            validity(t2)=th(t2).method(m).validity;
        end
        
        row=1;
        hs(row)=subplot(ROW,COL,row);
        hp=plot(x,pp); hold on
        hp.Color = color;
        hp.LineStyle = line;
        hp.Marker = marker;
        hs(row).YScale='log';
        hs(row).YLim=[1e-5 1];
        hs(row).YTick=[1e-4 1e-2 1];
        if m==1; hs(row).YLabel.String='p (t-test)'; end
        
        if m==size(th(t).method,2)
            xVal=4:0.5:7;
            hs(row).XTick = 1:length(xVal);
            hs(row).XTickLabel = {xVal};
            hs(row).XLim=[0.5 length(xVal)+0.5];
            hs(row).XLabel.String='Spike detection threshold factor';
            
            hp=plot([hs(row).XLim(1),hs(row).XLim(2)],[0.05 0.05],'--','Color',color1);
            hp=plot([hs(row).XLim(1),hs(row).XLim(2)],[0.01 0.01],'--','Color',color2);
            
            text(hs(row).XLim(2),0.05,'5% (*)','FontSize',8,'Color',color1)
            text(hs(row).XLim(2),0.01,'1% (**)','FontSize',8,'Color',color2)
            
            
            
            hs(row).YTick=[10^-4,10^-3,10^-2,10^-1,10^0];
            hs(row).YGrid = 'on';
            
            hl=legend(methodNames);
            hl.Location = 'Northoutside';
            
        end
        
    end
    
end

%% Nested Function

%% get plot styles
    function [color,line,marker]=getPlotStyle(idx)
        dark=.4;
        colors={[.4 .4 .4], [.4 .4 .4], [.4 .4 .4], [.6 .6 .6], [1 .1 0.5], [.1 .1 .5], [.8 .8 .8], [.8 .8 .8], [.8 .8 .8], [.8 .8 .8], [.9 .1 .1], [.9 .9 .9], [.8 .8 .8], [.7 .7 .7], [.6 .6 .6], [.5 .5 .5], [.4 .4 .4], [.3 .3 .3], [.2 .2 .2], [.1 .1 .1]};
        color=colors{idx};
        
        lines={'--', '-.', ':', '-', '--', '-.', ':', '-', '--', '--', '-.', ':', '-', '--', '-.', ':', '-', '--'};
        line= lines{idx};
        
        markers={'p', 'x', '*', 's', '.', 'd', '>', '^', '<', 'v','.','o', '^', 's', 'x', '<', '>', 'd', '*', 'p', 'v','.','o'};
        marker=markers{idx};
    end

end