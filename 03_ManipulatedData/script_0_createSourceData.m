%% Script: choose experimental data (and merge the 60s files to 300s files) that will be used in the next script as source for the in silico manipulated data
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
path_data = strrep(path_data,'03_ManipulatedData','01_CalcSync');
path_save = [path filesep 'Data_Source'];

list1=dir(path_data);

%% only for 10 �M BIC concentration (=last concentration)
for idl1=size(list1,1)
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
        MERGED_REF=Init_MERGED(); % needed for MergeTS function
        MERGED_BIC=Init_MERGED(); % needed for MergeTS function
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
            if f<=5 % ref (0 �M BIC)
                MERGED_REF=MergeTS(temp.SPIKEZ,MERGED_REF,fname,flag_realClock);
            else % 10 �M BIC
                MERGED_BIC=MergeTS(temp.SPIKEZ,MERGED_BIC,fname,flag_realClock);
            end
        end
        
        %% Save Results
        if ~exist(path_save,'dir')
            mkdir(path_save);
        end
        
        % Save whole chip
        TS = MERGED_REF.TS;
        TS(TS==0)=NaN;
        save([path_save filesep con(e).chip(c).name '_'  fname(1:end-4) '_REF.mat'],'TS')
        
        TS = MERGED_BIC.TS;
        TS(TS==0)=NaN;
        save([path_save filesep con(e).chip(c).name '_'  fname(1:end-4) '_BIC.mat'],'TS')
        
    end
end
disp('finished')