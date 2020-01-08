%% Script: Calculate synchrony values for all in silico manipulated data
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

%% path to data
dir_name=[path filesep 'Data_Shuffled'];

%% save path
folder_result = 'Result_Sync';

%% Synchrony measure methods
methods=    {   ...
                    'Sync_CC_selinger', ...
                    'Sync_MI1', ...
                    'Sync_MI2', ...
                    'Sync_PS', ...
                    'Sync_PS_M', ...
                    'Sync_Contrast_1ms', ...
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
method=methods([1,3,5,6,10,11,12,14,15]); % only relevant for paper


%% get all subdirectories
[dirarray,files]=subdir(dir_name);
if ~iscell(dirarray)
    dirarray = {dir_name}; % force it to be a cell array of strings in case only one file is selected
    tmp=dir(dir_name);
    files{1}={tmp.name};
end 


% for all subdirectories
for jjj=1:size(dirarray,2) 
    current_dir=dirarray{jjj};
    filearray=files{jjj};

    % for all files inside current directory
    for iii=1:size(filearray,2) % Loop from 1 to last selected file    
        current_file = filearray{iii}; 
        [~,filename,ext] = fileparts(current_file); % get file extension 
        if strcmp(ext,'.mat')
            
            % for all sync measure methods
            for m=1:size(method,2)
                
                % create new path and filenames
                idx = max(strfind(current_dir,[filesep 'Data_Shuffled' filesep])); % in case path has more than one '\Data\' in it, just change lowest hirachy 
                tmp_path = strrep(current_dir(idx:end),[filesep 'Data_Shuffled' filesep],[filesep folder_result filesep method{m} filesep]) ; % replace string: MODIFIEDSTR = strrep(ORIGSTR,OLDSUBSTR,NEWSUBSTR)
                folder_name = [current_dir(1:idx-1) tmp_path filesep]; % /Results_Sync/method/sourceFiles/shufflingMethod/n/percentage 
                filename=current_file;
                [~,fname,ext] = fileparts([folder_name filename]);
                fname_nu = [fname '.mat'];
               
                if ~exist([folder_name fname_nu],'file') % only calc if file not already exists
                   
                    display(method{m})                   
                    
                    %% Calculate synchrony 
                    S_tmp = load([current_dir filesep current_file]); % loads M_TS and percentage
                    S_tmp = S_tmp.S;
                    TS=S_tmp.M_TS;
                    
                    TS(TS==0)=NaN;
                    TS=sort(TS);

                    
                    % manipulated TS:
                    SPIKEZ.TS=TS;
                    SPIKEZ.AMP=zeros(size(TS));
                    SPIKEZ.PREF.SaRa = 10000;
                    SPIKEZ.PREF.rec_dur = 300;
                    time_win = 300;
                    FRmin=5;
                    [WIN,~]=CalcParameter_function(SPIKEZ,method{m},time_win,FRmin); % calculate all parameter
                    PARAMETER=unpackWIN2PARAMETER(WIN); % unpack structrue WIN and calculate electrode wise parameter 
                    S.SYNC = PARAMETER.mean;
                    S.percentage=S_tmp.percentage;
                    S.percentage_index = S_tmp.percentage_index;
                    
                    % original vs. random (same number of spikes!) ->
                    % needed for normalization of each method
                    TSr = TS;    
                    for el=1:size(TSr,2)
                        TSr(1:sum(~isnan(TS(:,el))),el)=rand(sum(~isnan(TS(:,el))),1) * time_win;
                    end
                    TSr=sort(TSr);
                    SPIKEZ.TS=TSr;
                    [WIN,~]=CalcParameter_function(SPIKEZ,method{m},time_win,FRmin); % calculate all parameter
                    PARAMETER=unpackWIN2PARAMETER(WIN); % unpack structrue WIN and calculate electrode wise parameter 
                    S.SYNC_originalVSrandom = PARAMETER.mean;
                    

                    %% save results
                    if ~exist(folder_name,'dir')
                        mkdir(folder_name);  
                    end
                    save([folder_name filename], 'S')
                end
            end
        end
    end
end

disp('finished')