%% Script: plot inter spike interval histogram of the experimental data (without BIC and with 10 µM BIC)  
% Paper title: "Comparison of different spike train synchrony measures regarding their robustness to erroneous data from bicuculline induced epileptiform activity"
% Author: Manuel Ciba (2019)
% 
% IMPORTANT: This scripts uses functions from "DrCell" (also included in
% this package). In order to use the DrCell functions just run the m-file DrCell.m. 
% This will add all DrCell functions to your Matlab path.

clear all
close all
clc

path_full=mfilename('fullpath'); % get path of this script
[path,~] = fileparts(path_full); % separate path from filename
cd(path)


%% Database path
path_data = [path filesep 'MEA-TS-Data' filesep 'MC_Th_Factor_5' filesep 'Data'];
path_pics = [path filesep 'MEA-TS-Data' filesep 'MC_Th_Factor_5' filesep 'Pics' filesep 'ISI'];



list1=dir(path_data);

%% for each BIC-Concentration (=experiment)
for idl1=3:size(list1,1)
    disp(['Concentration: ' list1(idl1).name])
    e=idl1-2; % index experiment
    con(e).name=list1(idl1).name;
    path2 = [path_data filesep list1(idl1).name];
    list2 = dir(path2);
    
    %% for each Chip
    for idl2=3:size(list2,1)
        disp(['   Chip: ' list2(idl2).name])
        c=idl2-2; % index chip
        con(e).chip(c).name=list2(idl2).name;
        path3 = [path2 filesep list2(idl2).name];
        list3 = dir(path3);
        
        
        %% load first 5 (=ref) and second 5 (=bic) files
        f=0; % index file
        MERGED=Init_MERGED(); % needed for MergeTS function
        for idl3=3:size(list3,1)
            
            [~,fname,ext] = fileparts(list3(idl3).name);
            
            % load TS.mat files:
            if strcmp(ext,'.mat')
                disp(['      File: ' list3(idl3).name])
                f=f+1; % index file
                temp_nu = load([path3 filesep list3(idl3).name]);    % load file
                temp.SPIKEZ=temp_nu.temp.SPIKEZ;
                temp.SPIKEZ.TS=temp.SPIKEZ.neg.TS;                   % here, TS are stored in field "neg"
                temp.SPIKEZ.AMP=temp.SPIKEZ.neg.AMP;                 % here, AMP are stored in field "neg"
            end
            
            % load Time-Amp txt files:
            if strcmp(ext,'.txt')
                disp(['      File: ' list3(idl3).name])
                f=f+1;
                data_raw = load([path3 filesep list3(idl3).name]);
                temp.SPIKEZ.TS = data_raw(2:end,1:2:end);
                temp.SPIKEZ.AMP = data_raw(2:end,2:2:end);
                temp.SPIKEZ.PREF.rec_dur=60;
                temp.SPIKEZ.PREF.SaRa=10000;
                temp.SPIKEZ.AMP(temp.SPIKEZ.TS==0)=NaN; % NaN-padding
                temp.SPIKEZ.TS(temp.SPIKEZ.TS==0)=NaN; % NaN-padding
            end
            
            con(e).chip(c).file(f).name=list3(idl3).name;
            con(e).chip(c).file(f).SPIKEZ=temp.SPIKEZ;
            
            % Merge all 10 fiels together (=10 minutes) and seperate them
            % later into two 5 minute files (=first 5 minutes: ref, rest:
            % BIC)
            flag_realClock=0; % just combine single TS files by merging it toghether without 1 second break after each recording
            fname=list3(idl3).name;
            MERGED=MergeTS(temp.SPIKEZ,MERGED,fname,flag_realClock);
        end
        con(e).chip(c).MERGED=MERGED;
    end
end


%% plot
if 0     % plot IBI vs BD vs SIB for every chip
    vm = 10;%[0.1 0.2 0.5 1 2 5 10]; % bic concentrations
    Rows = 2;%size(vm,2);
    allChips = {'19498','20692','24370','24659','25165'}; %{'18812','19498','19579','20692','24370','24659','25050','25165'};
    Cols=size(allChips,2); % number of chips
    
    row = 1;
    col=1;
    
    hf=figure('Name','1','visible','on');
    
    for e=size(con,2) % only for 10 uM
        
        for c=1:size(con(e).chip,2) % for every Chip
            
            [~,col]=max(strcmp(allChips,con(e).chip(c).name));
            row=1;
            
            % Burstdetection to get IBI and BD
            Selection={ ...
                'BR_baker100',...
                'BD_baker100',...
                'SIB_baker100',...
                'IBI_baker100',...
                'BR_selinger',...
                'BD_selinger',...
                'SIB_selinger',...
                'IBI_selinger'
                };
            time_win = 300; % here, the 600s file is splitted into two (ref and BIC)
            FR_min = 5; % at least 5 spikes have to be on electrode in order to consider it active
            [WIN,MERGED]=CalcParameter_function(con(e).chip(c).MERGED,Selection,time_win,FR_min);
            
            for el=1:60
                REF(c).BD(el)=median(nonzeros(WIN(1).parameter(6).values(:,el))); % BD Selinger
                REF(c).SIB(el)=median(nonzeros(WIN(1).parameter(7).values(:,el))); % SIB Selinger
                REF(c).IBI(el)=median(nonzeros(WIN(1).parameter(8).values(:,el))); % IBI Selinger
                
                BIC(c).BD(el)=median(nonzeros(WIN(2).parameter(6).values(:,el))); % BD Selinger
                BIC(c).SIB(el)=median(nonzeros(WIN(2).parameter(7).values(:,el))); % SIB Selinger
                BIC(c).IBI(el)=median(nonzeros(WIN(2).parameter(8).values(:,el))); % IBI Selinger
            end
            
            % Plot:
            
            
            % Ref:
            row=1;
            hs = subplot(Rows, Cols, Cols*(row-1)+col);
            hp=stem3(REF(c).IBI, REF(c).BD, REF(c).SIB); hold on
            % Make plot quadradic
            maxLim = max([hp.Parent.XLim hp.Parent.YLim]);
            maxLim = 15;
            hp.Parent.XLim = [0 maxLim];
            hp.Parent.YLim = [0 maxLim];
            plot([0 maxLim],[0 maxLim],'k--'); % plot diagonal
            hs.XLabel.String='IBI /s';
            hs.YLabel.String='BD /s';
            hs.ZLabel.String='SIB';
            hs.Title.String=con(e).chip(c).name;
            
            % BIC:
            row=2;
            hs = subplot(Rows, Cols, Cols*(row-1)+col);
            hp=stem3(BIC(c).IBI, BIC(c).BD, BIC(c).SIB); hold on
            % Make plot quadradic
            maxLim = max([hp.Parent.XLim hp.Parent.YLim]);
            maxLim = 15;
            hp.Parent.XLim = [0 maxLim];
            hp.Parent.YLim = [0 maxLim];
            plot([0 maxLim],[0 maxLim],'k--'); % plot diagonal
            hs.XLabel.String='IBI /s';
            hs.YLabel.String='BD /s';
            hs.ZLabel.String='SIB';
            
            
            
            
            %plotISI_Histogram(getLogISI(con(e).chip(c).MERGED.TS));
            %hs=plotSpikeTrain(con(e).chip(c).MERGED.TS, hs, 'dot', 1);
            %hs.XLabel.String=[num2str(vm(1)) 'uM BIC'];
            
            
        end
        
    end
    
    if 1
        if ~exist(path_pics,'dir')
            mkdir(path_pics);
        end
        
        hf.PaperPositionMode = 'auto';
        hf.Units='points';
        hf.Position=[0,0,1920,1080];
        print([path_pics filesep 'ISI'],'-dpng','-r600')
    end
    
end

%% plot
if 1     % plot IBI vs BD vs SIB from all chips in one plot AND ISI-Histogram
    vm = 10;%[0.1 0.2 0.5 1 2 5 10]; % bic concentrations
    Rows = 2;%size(vm,2);
    allChips = {'19498','20692','24370','24659','25165'}; %{'18812','19498','19579','20692','24370','24659','25050','25165'};
    Cols=size(allChips,2); % number of chips
    
    row = 1;
    col=1;
    
    hf=figure('Name','1','visible','on');
    
    % init variables
    REF.logISI = 0;
    BIC.logISI = 0;
    
    for e=size(con,2) % only for 10 uM
        
        for c=1:size(con(e).chip,2) % for every Chip
            
            [~,col]=max(strcmp(allChips,con(e).chip(c).name));
            row=1;
            
            % Burstdetection to get IBI and BD
            Selection={ ...
                'BR_baker100',...
                'BD_baker100',...
                'SIB_baker100',...
                'IBI_baker100',...
                'BR_baker200',...
                'BD_baker200',...
                'SIB_baker200',...
                'IBI_baker200'
                };
            time_win = 300; % here, the 600s file is splitted into two (ref and BIC)
            FR_min = 5; % at least 5 spikes have to be on electrode in order to consider it active
            [WIN,MERGED]=CalcParameter_function(con(e).chip(c).MERGED,Selection,time_win,FR_min);
            
            for el=1:60
                BR=WIN(1).parameter(5).allEl(el); % BR Selinger per minute
                if BR >= 1 % if more than one burst per minute consider el. active
                    REF.BD(c,el)=median(nonzeros(WIN(1).parameter(6).values(:,el))); % BD Selinger
                    REF.SIB(c,el)=median(nonzeros(WIN(1).parameter(7).values(:,el))); % SIB Selinger
                    REF.IBI(c,el)=median(nonzeros(WIN(1).parameter(8).values(:,el))); % IBI Selinger
                else
                    REF.BD(c,el)=NaN; % BD Selinger
                    REF.SIB(c,el)=NaN; % SIB Selinger
                    REF.IBI(c,el)=NaN; % IBI Selinger
                end
                
                BR=WIN(2).parameter(5).allEl(el); % BR Selinger per minute
                if BR >= 1 % if more than one burst per minute consider el. active
                    BIC.BD(c,el)=median(nonzeros(WIN(2).parameter(6).values(:,el))); % BD Selinger
                    BIC.SIB(c,el)=median(nonzeros(WIN(2).parameter(7).values(:,el))); % SIB Selinger
                    BIC.IBI(c,el)=median(nonzeros(WIN(2).parameter(8).values(:,el))); % IBI Selinger
                else
                    BIC.BD(c,el)=NaN; % BD Selinger
                    BIC.SIB(c,el)=NaN; % SIB Selinger
                    BIC.IBI(c,el)=NaN; % IBI Selinger
                end
            end
            
            
            % safe ISIs of all data
            logISI = getLogISI(WIN(1).SPIKEZ.TS);
            REF.logISI(end+1:end+length(logISI)) = logISI;
            logISI = getLogISI(WIN(2).SPIKEZ.TS);
            BIC.logISI(end+1:end+length(logISI)) = logISI;
            
        end
        
        
        % Plot ISI-histogram
        hs = subplot(2,1,1);
        hp = plotISI_Histogram(REF.logISI,hs); hold on
        hp.FaceColor = 'k';
        hp.EdgeAlpha = 0; % make edges invisible
        hp.FaceAlpha = 0.5; % make color transparent
        hp = plotISI_Histogram(BIC.logISI,hs); hold on
        hp.FaceColor = 'g';
        hp.EdgeAlpha = 0; % make edges invisible
        hp.FaceAlpha = 0.5; % make color transparent
        legend('0 �M BIC','10 �M BIC')
       
        % Plot (IBI vs BD vs SIB):
        hs = subplot(2, 2, 3);
        hp1=stem3(REF.IBI(:),REF.BD(:),REF.SIB(:),'k'); hold on
        hp2=stem3(BIC.IBI(:),BIC.BD(:),BIC.SIB(:),'g'); hold on
        hs.XLabel.String='IBI /s';
        hs.YLabel.String='BD /s';
        hs.ZLabel.String='SIB';
        maxLim = 15;
        hs.XLim = [0 maxLim];
        hs.YLim = [0 maxLim];
        plot([0 maxLim],[0 maxLim],'k--'); % plot diagonal
        legend([hp1, hp2],'0 �M BIC','10 �M BIC')
        
        % Plot (IBI/BD vs SIB)
        hs = subplot(2,2,4);
        REFratio = REF.IBI(:)./REF.BD(:);
        BICratio = BIC.IBI(:)./BIC.BD(:);
        hp=plot(REFratio,REF.SIB(:),'kx'); hold on
        hp=plot(BICratio,BIC.SIB(:),'go'); hold on
        hs.XLabel.String='IBI/BD';
        hs.YLabel.String='SIB';
        hs.XScale='log';
        legend('0 �M BIC','10 �M BIC')
        
    end
    
    if 1
        if ~exist(path_pics,'dir')
            mkdir(path_pics);
        end
        
        hf.PaperPositionMode = 'auto';
        hf.Units='centimeter';
        hf.Position=[0,0,20,15];
        print([path_pics filesep 'IBIvsBDvsSIB_allInOne'],'-dpng','-r600')
    end
    
end





disp('finished')
