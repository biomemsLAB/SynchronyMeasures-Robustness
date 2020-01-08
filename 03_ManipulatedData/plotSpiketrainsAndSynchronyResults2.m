function plotSpiketrainsAndSynchronyResults2(Struc)

% Titles of manipulation methods
Titles={{'1) Added spikes' '(Requirement: Robust)'},{'2) Deleted spikes' '(Requirement: Robust)'}};

% Size of Figure
Units = 'centimeter';
Hi = 23*1.3 *(5/6);
W = 18*1.3;

exampleSpikeTrainNumber = 4;
flag_showCurves = 0; % 1: show synchrony over manipultion level curve, 0: don't show it

for s=1:size(Struc.SyncMethod,2)
    disp(Struc.SyncMethod(s).name)
    for f=1:1%size(Struc.SyncMethod(s).Testfile,2)
        
        % only create figures at first time -> otherwise it will overwrite it
        if s==1
            hf(f) = figure(f);
            hf(f).Units = Units;
            hf(f).Position = [0,0,W,Hi]; % change size of displayed window
            hf(f).PaperPosition=[0,0,W,Hi]; % change size of file when "print" is used
        end
        
        for m=1:size(Struc.SyncMethod(s).Testfile(f).ShufflingMethod,2)
            %% unpack structure:
            percentage = Struc.percentage;
            TS = Struc.SyncMethod(s).ShufflingMethod(m).Testfile(exampleSpikeTrainNumber).TS(1).TS;
            Sync_raw = Struc.SyncMethod(s).ShufflingMethod(m).Sync_raw; % f x n x p
            Sync_originalVSrandom = Struc.SyncMethod(s).ShufflingMethod(m).Sync_originalVSrandom;
            
            % normalize all synchrony values for each chip to its first
            % value (so every synchrony values starts with 1). Dimensions:
            % f x n x p (=chip x number of simulations x level of manipulation)
            Sync_raw_scaled = zeros(size(Sync_raw));
            Sync_raw_norm = zeros(size(Sync_raw));
            for iii = 1:size(Sync_raw,1) % for every chip
                for p=1:size(Sync_raw,3) % for every percentage
                    
                    % 1) range-correct as every sync measures has different
                    % range -> makes it possible to compare sync measures.
                    % scale Sync values between 0 and 1: (Sync-min)/(max-min) = (Sync-random)/(1-random)
                    Sync_raw_scaled(iii,:,p) = (Sync_raw(iii,:,p)-mean(Sync_originalVSrandom(iii,:,p)))./(1-mean(Sync_originalVSrandom(iii,:,p)));
                    
                    % 2) normalize to first value (unmanipulated data) as every chip/network shows different synchrony levels.
                    % -> makes it possible to compare different chips
                    % (=calculate mean over all chips)
                    Sync_raw_norm(iii,:,p) = Sync_raw_scaled(iii,:,p)./mean(Sync_raw_scaled(iii,:,1));
                end
            end
            
            % rearrange f x n x p to        f and n x p
            for p=1:size(Sync_raw,3)
                tmp = Sync_raw(:,:,p);
                tmp_Sync_raw(:,p) = tmp(:);
                tmp = Sync_raw_scaled(:,:,p);
                tmp_Sync_raw_scaled(:,p) = tmp(:);
                
                
                tmp =  Sync_originalVSrandom(:,:,p);
                tmp_originalVSrandom(:,p)= tmp(:);
                tmp = Sync_raw_norm(:,:,p);
                tmp_Sync_raw_norm(:,p) = tmp(:);
            end
            Sync_raw = tmp_Sync_raw;
            Sync_raw_scaled = tmp_Sync_raw_scaled;
            Sync_originalVSrandom = tmp_originalVSrandom;
            Sync_raw_norm = tmp_Sync_raw_norm;
            
            %% calculate mean std ect. of sync values
            Sync_mean = mean(Sync_raw,'omitnan'); % n x percentage
            Sync_std = std(Sync_raw,'omitnan');
            Sync_min = min(Sync_raw,[],'omitnan');
            Sync_max = max(Sync_raw,[],'omitnan');
            Sync_median = median(Sync_raw,'omitnan');
            
            Sync_mean_norm = mean(Sync_raw_norm,'omitnan'); % n x percentage
            Sync_std_norm = std(Sync_raw_norm,'omitnan');
            Sync_min_norm = min(Sync_raw_norm,[],'omitnan');
            Sync_max_norm = max(Sync_raw_norm,[],'omitnan');
            Sync_median_norm = median(Sync_raw_norm,'omitnan');
            
            MethodRank(s).Name = Struc.SyncMethod(s).name;
            MethodRank(s).sum(m) = sum(std(Sync_raw_norm));
            
            %% plot
            hf(f)=figure(hf(f)); % swtich between different figures -> make figure hf(f) the current figure
            col = m;
            
            % REF Spiketrain
            row=1;
            if s==1
                if flag_showCurves
                    ROW = 6; % 1st row full ref spiketrain, 2nd row TS, 3rd row Sync, 4th row scaled. synchrony, 5th empty place for legend, 6th summary of results
                else
                    ROW = 3;
                end
                hs(row,1) = subplot(ROW,1,1);
                TS = Struc.SyncMethod(1).Testfile(exampleSpikeTrainNumber).ShufflingMethod(col).TS(1).TS;
                hs(row,1) = plotSpikeTrain(TS,hs(row,1),'dot');
                hs(row,1).XLim=[0 15];
                hs(row,1).YLabel.String={'Spike trains' '(Manipulation' 'level: 0)'};
                hs(row,1).Title.String={'Recorded spike trains'};
            end
            
            COL=size(Struc.SyncMethod(s).Testfile(f).ShufflingMethod,2); % shuffling methods
            
            % TS
            row=row+1;
            if s==1 % only plot TS one time
                TS = Struc.SyncMethod(1).Testfile(exampleSpikeTrainNumber).ShufflingMethod(col).TS(11).TS;
                hs(row,col) = subplot(ROW,COL,col+2);
                hs(row,col) = plotSpikeTrain(TS,hs(row,col),'dot',0);
                hs(row,col).XLim = [0 15];
                string1 = '(Manipulation';
                string2 = ['level: ' , num2str(percentage(11)), ')'];
                hs(row,col).YLabel.String={'Spike trains', string1, string2};
                hs(row,col).Title.String = Titles{m};
                
                % Change size and position of first plot (=ref spike trains)
                hs(1,1).Position = [0.5-hs(2,1).Position(3)/2 hs(1,1).Position(2)*1.01 hs(2,1).Position(3) hs(1,1).Position(4)];
            end
            
            % Sync
            if flag_showCurves
                row = row+1;
                hs(row,col) = subplot(ROW,COL,COL*(row-1)+col);
                x=percentage;
                y=Sync_mean;
                dmin=Sync_std;
                dmax=Sync_std;
                [color]=getPlotStyle(s);
                hold all;
                h(row,col)=errorbar(x,y,dmin,dmax);
                [h(row,col).Color,h(row,col).LineStyle,h(row,col).Marker]=getPlotStyle(s);
                hs(row,col).XLim=[0 1];
                hs(row,col).XScale='lin';
                hs(row,col).YScale='lin';
                hs(row,col).YLim=[0 1];
                if col==1; hs(row,col).YLabel.String = {'Absolute' 'synchrony'}; end
            end
            
            % Sync (normalized)
            if flag_showCurves
                row = row+1;
                hs(row,col) = subplot(ROW,COL,COL*(row-1)+col); 
                x=percentage;
                y=Sync_mean_norm;
                dmin=Sync_std_norm;
                dmax=Sync_std_norm;
                [color]=getPlotStyle(s);
                hold all;
                h(row,col)=errorbar(x,y,dmin,dmax);
                %h(row,col)=plot(x,y);
                [h(row,col).Color,h(row,col).LineStyle,h(row,col).Marker]=getPlotStyle(s);
                hs(row,col).XLim=[0 1];
                hs(row,col).XScale='lin';
                hs(row,col).YLim=[-inf inf];
                hs(row,col).XLabel.String = 'Manipulation level';
                if col==1; hs(row,col).YLabel.String = {'Normalized' 'synchrony'}; end
            end
            
            
        end % end ShufflingMethod
    end % end Testfile
    
    MethodRank(s).finalSum = sum(MethodRank(s).sum);
end % end SyncMethod



%% Legend
customNames = {'CC (bin = 500 ms)', 'MI (bin = 500 ms)', 'STTC (dt = 100 ms)', 'PS', 'Spike-contrast', 'A-SPIKE-synchronization', 'A-ISI-distance', 'A-SPIKE-distance', 'ARI-SPIKE-distance'};
if flag_showCurves
    for i = 1:length(Struc.SyncMethod)
        Struc.SyncMethod(i).name = strrep(Struc.SyncMethod(i).name,'Sync_',''); 
        Struc.SyncMethod(i).name = strrep(Struc.SyncMethod(i).name,'_',''); 
    end
    axes(hs(4,2))
    hl=legend(customNames,'location','north');
    hl.Position=[.5 - hl.Position(3)/2, hs(4,1).Position(2) - hl.Position(4)*1.2 , hl.Position(3), hl.Position(4)];
    hl.AutoUpdate = 'off';
    pause(2)
    hs(4,2).Position=[hs(4,2).Position(1), hs(4,1).Position(2), hs(3,2).Position(3), hs(3,2).Position(4)]; % make legend subplot same size as other plots (size is changed due to legend)
end

%% lines of spike-contrast to top layer
if flag_showCurves
    hs(3,1).Children = hs(3,1).Children([5 1 2 3 4 6 7 8 9]);
    hs(4,1).Children = hs(4,1).Children([5 1 2 3 4 6 7 8 9]);  
    hs(3,2).Children = hs(3,2).Children([5 1 2 3 4 6 7 8 9]);
    hs(4,2).Children = hs(4,2).Children([5 1 2 3 4 6 7 8 9]);
end

%% Rank
for col=1:2
    row=ROW;
    
    if flag_showCurves
        hs(row,col)=subplot(ROW,COL,10+col);
    else
        hs(row,col)=subplot(ROW,COL,4+col);
    end
    
    
    for i=1:size(MethodRank,2)
        percentageValue(i) = MethodRank(i).sum(col);
        robustness(i) = -1 * percentageValue(i);
    end
    
    [sortedValues, index] = sort(percentageValue);
    
    ii=0;
    for i=index
        ii=ii+1;
        [color,line,marker]=getPlotStyle(i);
        
        hp=bar(ii,percentageValue(i)); hold all
        hp.LineStyle = line;
        hp.FaceColor = color;
    end
    hs(row,col).XLabel.String = [];
    if col==1
        hs(row,col).YLabel.String = {'Total deviation of' 'normalized synchrony'};
    else
        hs(row,col).YLabel.String = {''};
    end
    hs(row,col).XTick = [1:length(customNames)];
    hs(row,col).XTickLabel = customNames(index);
    hs(row,col).XTickLabelRotation = 45;
    hs(row,col).YLim=[0 20];
    box off
end

%% (a) (b) (c) ...
if flag_showCurves
    array = [1 2 3 4 6];
else
    array = 1:ROW;
end
for i= array
    x = hs(i,1).Position(1);
    y = hs(i,1).Position(2);
    h = hs(i,1).Position(4);
    
    if i==6; i=5; end
    
    ha=annotation('textbox', [x-0.11, y+h+0.01, 1, 0], 'string', ['(',char(i+96),')'],'color','k','fontw','b');
    ha.FontSize=10;
    ha.LineStyle='none';
end



%% Nested Function

%% get plot styles
    function [color,line,marker]=getPlotStyle(idx)
        dark=.4;
        % green colors={[.4 .4 .4], [.4 .4 .4], [.4 .4 .4], [.6 .6 .6], [.0 .8 0.2], [.8 .8 .8], [.8 .8 .8], [.8 .8 .8], [.8 .8 .8], [.8 .8 .8], [.9 .1 .1], [.9 .9 .9], [.8 .8 .8], [.7 .7 .7], [.6 .6 .6], [.5 .5 .5], [.4 .4 .4], [.3 .3 .3], [.2 .2 .2], [.1 .1 .1]};
        colors={[.4 .4 .4], [.4 .4 .4], [.4 .4 .4], [.6 .6 .6], [1 .1 0.5], [.1 .1 .5], [.8 .8 .8], [.8 .8 .8], [.8 .8 .8], [.8 .8 .8], [.9 .1 .1], [.9 .9 .9], [.8 .8 .8], [.7 .7 .7], [.6 .6 .6], [.5 .5 .5], [.4 .4 .4], [.3 .3 .3], [.2 .2 .2], [.1 .1 .1]};
        color=colors{idx};
        
        lines={'--', '-.', ':', '-', '--', '-.', ':', '-', '--', '--', '-.', ':', '-', '--', '-.', ':', '-', '--'};
        line= lines{idx};
        
        markers={'p', 'x', '*', 's', '.', 'd', '>', '^', '<', 'v','.','o', '^', 's', 'x', '<', '>', 'd', '*', 'p', 'v','.','o'};
        marker=markers{idx};
    end

end