%% Script: Calculate synchrony using different synchrony measueres to analyze the experimental data (without BIC and with 10 µM BIC)  
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
path_results = [path filesep 'MEA-TS-Data' filesep 'MC_Th_Factor_5' filesep 'Results'];


%% Select Parameter 
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
                    
                    
 methods=AllParameter([18:end]); % only apply synchrony measures


%% for each Method
for idl0=1:size(methods,2)
    disp(methods{idl0})
    pause(0);
    list1=dir(path_data);
    
    %% for each BIC-Concentration (=experiment)
    for idl1=3:size(list1,1)
       disp(['Concentration: ' list1(idl1).name])
       e=idl1; % index experiment
       exp(e).name=list1(idl1).name;
       path2 = [path_data filesep list1(idl1).name];
       list2 = dir(path2);

       %% for each Chip
       for idl2=3:size(list2,1)
           disp(['   Chip: ' list2(idl2).name])
           c=idl2-2; % index chip
           exp(e).chip(c).name=list2(idl2).name;
           path3 = [path2 filesep list2(idl2).name];
           list3 = dir(path3);
           
           folder_name = [path_results filesep methods{idl0} filesep list1(idl1).name filesep list2(idl2).name]; % method/concentration/chip/ 

           
           if ~exist([folder_name filesep 'Parameter.mat'],'file') % only calc if file not already exists
               
               %% for all files (1-5: ctrl, 6-10: bic, 11-15: washout (not available for all chips))
               f=0; % index file
               MERGED=Init_MERGED(); % needed for MergeTS function
               for idl3=3:12%size(list3,1)

                   [~,fname,ext] = fileparts(list3(idl3).name);
                   fname_nu = [fname '.mat'];
                       %% load TS files

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
                           temp.SPIKEZ.AMP(temp.SPIKEZ.TS==0)=NaN; % NaN-patting
                           temp.SPIKEZ.TS(temp.SPIKEZ.TS==0)=NaN; % NaN-patting
                       end

                       flag_realClock=0; % just combine single TS files by merging it toghether without 1 second break after each recording
                       MERGED=MergeTS(temp.SPIKEZ,MERGED,fname,flag_realClock); 
               end
           
               %% calculate all parameter
               FRmin=6; % only electrodes with more or equal than 6 spikes per minute (=min. firing rate=FRmin) are considered
               time_win= 300; % here for 300 seconds = 5 minutes 
               MERGED.TS(isnan(MERGED.TS))=0; % ZERO-PATTING!!!
               [WIN,MERGED]=CalcParameter_function(MERGED,methods{idl0},time_win,FRmin); % calculate all parameter
               PARAMETER=unpackWIN2PARAMETER(WIN); % unpack structrue WIN and calculate electrode wise parameter 
               PARAMETER=rmfield(PARAMETER,'values'); % remove values as they contain complete TS file in case for parameter "spikerate"

               %% Save Results
               % create folder for sync. measure results
               if ~exist(folder_name,'dir')
                mkdir(folder_name);  
               end
               save([folder_name filesep 'Parameter.mat'],'PARAMETER')
           
           end % end: if ~exist(parameter.mat)
           
       end
    end
end
disp('finished')
