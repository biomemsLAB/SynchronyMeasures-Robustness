%% Script: Analyze and plot results of the calculated synchrony values of the experimental data with different spike detection threshold factors
% Paper title: "Comparison of different spike train synchrony measures regarding their robustness to erroneous data from bicuculline induced epileptiform activity"
% Author: Manuel Ciba (2019)
%
% IMPORTANT: This scripts uses functions from "DrCell" (also included in
% this package). In order to use the DrCell functions just run the m-file DrCell.m.
% This will add all DrCell functions to your Matlab path.

%% CHANGE NAME OF FOLDER "100" to "99" so order of plot will be correct

clear all
close all
clc

path_full=mfilename('fullpath'); % get path of this script
[path,~] = fileparts(path_full); % separate path from filename
cd(path)

%% Database path
path_results = [path filesep 'MEA-TS-Data' filesep 'MC_Th_Factor_4' filesep 'Results_VarThresholds'];


%% Electrode wise parameter or normal parameter?
flag_electrodewise = 0;
if flag_electrodewise 
    P_index=2;
else
    P_index=1;
end

%% Select SyncMeasure Methods
AllParameter={ ...
                    'Spikerate', ...
                    'Amplitude', ...
                    'ActiveElectrodes', ...
                    'BR_baker100', ...
                    'BD_baker100', ...
                    'SIB_baker100', ...
                    'IBI_baker100', ...
                    'BR_selinger', ...
                    'BD_selinger', ...
                    'SIB_selinger', ...
                    'IBI_selinger', ...
                    'NBR_chiappalone', ...
                    'NBD_chiappalone', ...
                    'SINB_chiappalone', ...
                    'INBI_chiappalone', ...
                    'Entropy_bin100', ...
                    'Entropy_capurro', ...
                    'Sync_CC_selinger', ...
                    'Sync_MI1', ...
                    'Sync_MI2', ...
                    'Sync_PS', ...
                    'Sync_PS_M', ...
                    'Sync_Contrast', ...
                    'Sync_ISIDistance', ...
                    'Sync_SpikeDistance', ...
                    'Sync_SpikeSynchronization', ...
                    'Sync_ASpikeSynchronization', ...
                    'Sync_AISIDistance', ...
                    'Sync_ASpikeDistance', ...
                    'Sync_RISpikeDistance', ...
                    'Sync_RIASpikeDistance', ...
                    'Sync_STTC'
                    };    

methods=AllParameter([18,20,32,22, 23, 27,28,29,31]); % only this methods are used for the paper
methodNames = { 'CC (bin = 500 ms)', 'MI (bin = 500 ms)', 'STTC (dt = 100 ms)', 'PS', 'Spike-contrast', 'A-SPIKE-synchronization', 'A-ISI-distance', 'A-SPIKE-distance', 'ARI-SPIKE-distance'};

%% 1) Load Sync-Data

list0 = dir(path_results);

for idl0=3:size(list0,1)
    
    disp(['Threshold: ' list0(idl0).name])
    t=idl0-2;
    th(t).name = list0(idl0).name;

    % for each Method
    for idl1=1:size(methods,2)
       disp(['Method: ' methods{idl1}])
       m=idl1; % index method
       th(t).method(m).name=methods{idl1};
       path2 = [path_results filesep list0(idl0).name filesep methods{idl1}];
       list2 = dir(path2);

       % for each Experiment
       for idl2=size(list2,1) %% only last BIC concentration!!!!
           disp(['   Experiment: ' list2(idl2).name])
           e=1; % index experiment
           th(t).method(m).exp(e).name=list2(idl2).name;
           path3 = [path2 filesep list2(idl2).name];
           list3 = dir(path3);

           % for each Chip
           for idl3=3:size(list3,1)
               disp(['      Chip: ' list3(idl3).name])
               c=idl3-2; % index method
               th(t).method(m).exp(e).chip(c).name=list3(idl3).name;          
               path4 = [path3 filesep list3(idl3).name];
               list4 = dir(path4);

               % Parameter.mat (first 300 s: ref (PARAMETER.mean(1)), second 300 s: bic (PARAMETER.mean(2)))
               temp=load([path4 filesep 'Parameter.mat']);
               th(t).method(m).exp(e).chip(c).PARAMETER = temp.PARAMETER;
               clear temp
           end
       end
    end
end

%% 2) Calculate mean values
for t=1:size(th,2)
    for m=1:size(th(t).method,2)
       for e=1

          REF_c=zeros(1,size(th(t).method(m).exp(e).chip,2));
          BIC_c=zeros(1,size(th(t).method(m).exp(e).chip,2));
          DIF_c=zeros(1,size(th(t).method(m).exp(e).chip,2));

          for c=1:size(th(t).method(m).exp(e).chip,2)

              % get mean parameter value for every chip
              REF_c(c) = th(t).method(m).exp(e).chip(c).PARAMETER(P_index).mean(1);
              BIC_c(c) = th(t).method(m).exp(e).chip(c).PARAMETER(P_index).mean(2);
              DIF_c(c) = BIC_c(c) / REF_c(c);

          end

          % calculate mean value of all chips
          th(t).method(m).REFmean(e) = mean(REF_c);
          th(t).method(m).BICmean(e) = mean(BIC_c);
          th(t).method(m).DIFmean(e) = mean(DIF_c);
          th(t).method(m).REF = REF_c;
          th(t).method(m).BIC = BIC_c;
          % calculate std value of all chips
          th(t).method(m).REFstd(e) = std(REF_c);
          th(t).method(m).BICstd(e) = std(BIC_c);
          th(t).method(m).DIFstd(e) = std(DIF_c);

          [pp,ppf,ppf_rel,validity]=pairedTest(BIC_c,REF_c,'right'); % Tail='both', test if BIC!=REF, Tail='right', test if BIC > REF
          th(t).method(m).pp(e)=pp;      % parameter dependend test
          th(t).method(m).validity(e)=validity;
          th(t).method(m).ppf(e)=ppf; % parameter free test
          th(t).method(m).ppf_rel(e)=ppf_rel; % parameter free test using relative values
       end
    end
end

%% 3) Plot Sync-Values
if 1
    plotSyncValuesVsBic_VarTh(th,methodNames)
end