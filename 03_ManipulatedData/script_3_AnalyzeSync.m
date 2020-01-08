%% Script: Analyze and plot the in silico manipulated data 
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
SyncMethods=methods([1, 3, 15, 5, 6, 10, 11, 12, 14 ]); % only this methods are used for the paper

if ~exist([path filesep 'temp_Struc.mat'],'file') % in temp_Struc.mat all calculated synchrony values are saved


    %% Preferences
    folder_TS = 'Data_Shuffled';
    folder_Sync = 'Result_Sync';
    path_TS = [path filesep folder_TS]; % Data_Shuffled/SourceFile/ShuffleMethod/n/percentage.mat
    path_Sync = [path filesep folder_Sync];
    N=1; % here only selected N is displayed (1...20)
    rec_dur=300;


    path0 = path_Sync;
    list0 = dir(path0); 

    %% for every Sync-Measure method
    s=0;
    for s=1:size(SyncMethods,2)
        path1 = [path0 filesep SyncMethods{s}];
        list1 = dir(path1);
        Struc.SyncMethod(s).name = SyncMethods{s};

        %% for every Test-File
        f=0;
        for idx1=3:size(list1,1) 
            disp(['Testfile: ' list1(idx1).name])
            f=f+1;
            path2 = [path1 filesep list1(idx1).name];
            list2 = dir(path2);
            Struc.SyncMethod(s).Testfile(f).name = list1(idx1).name;

            %% for every shuffling method
            m=0;
            for idx2=3:size(list2,1) 
                m=m+1;
                path3 = [path2 filesep list2(idx2).name];
                list3 = dir(path3);
                Struc.SyncMethod(s).Testfile(f).ShufflingMethod(m).name = list2(idx2).name;

                %% load TS
                % for every percentage value
                for idx3=3:size(list3,1)
                    p=idx3-2;
                    path4 = [path3 filesep list3(idx3).name];
                    path4_TS = strrep(path4,[folder_Sync filesep SyncMethods{s}],folder_TS); % Switch from Sync-Database to TS-Database, NOTE: sync-Database also includes Sync-Measures in path
                    list4 = dir(path4_TS);
                    clear S
                    temp = load([path4_TS filesep list4(2+N).name]); % load Nth sample
                    S=temp.S;
                    TS(S.percentage_index).TS=S.M_TS; 
                end % END: for all percentages x 

                %% load Sync Values
                % for every percentage value
                for idx3=3:size(list3,1)
                    p=idx3-2; % index percentage
                    path4 = [path3 filesep list3(idx3).name];
                    list4 = dir(path4);

                    % for all n
                    for idx4=3:size(list4,1) % 22 because in some folders there is a file called 60.mat which shouldn't be used
                        n=idx4-2; % index sample
                        temp = load([path4 filesep list4(idx4).name]); % load Nth sample
                        S=temp.S;
                        Sync_raw(n,S.percentage_index) = S.SYNC; % synchrony value from original vs. manipulated
                        Sync_originalVSrandom(n,S.percentage_index) = S.SYNC_originalVSrandom; % synchrony value from original vs. random signal
                        percentage(S.percentage_index)=S.percentage;
                    end % END: for all n  

                end % END: for all percentages x 

                %% fill structure:
                Struc.SyncMethod(s).Testfile(f).ShufflingMethod(m).TS = TS; % dim: n x percentage_manipulation
                Struc.SyncMethod(s).Testfile(f).ShufflingMethod(m).Sync_raw = Sync_raw;
                Struc.SyncMethod(s).Testfile(f).ShufflingMethod(m).Sync_originalVSrandom = Sync_originalVSrandom;
                Struc.percentage = percentage;

            end % END: for all shuffling methods
        end % END: for all Test files
        
        % rearrange structure so all testfiles are merged together like
        % this: Struc.SyncMethod(s).ShufflingMethod(m).Sync(f,n,percentage)
        for m=1:size(Struc.SyncMethod(s).Testfile(1).ShufflingMethod,2)
           for f=1:size(Struc.SyncMethod(s).Testfile,2)
               Struc.SyncMethod(s).ShufflingMethod(m).Testfile(f).TS=Struc.SyncMethod(s).Testfile(f).ShufflingMethod(m).TS;
               Struc.SyncMethod(s).ShufflingMethod(m).Sync_raw(f,:,:)=Struc.SyncMethod(s).Testfile(f).ShufflingMethod(m).Sync_raw;
               Struc.SyncMethod(s).ShufflingMethod(m).Sync_originalVSrandom(f,:,:)=Struc.SyncMethod(s).Testfile(f).ShufflingMethod(m).Sync_originalVSrandom;
           end
        end
        
        % delete old fields:
        %Struc=rmfield(Struc.SyncMethod(s),'Testfile');
        
    end % END: for all sync-measure 

    save('temp_Struc.mat','Struc','-v7.3') % save structure to avoid long loading times

else
    load('temp_Struc.mat')
end

%% Plot
plotSpiketrainsAndSynchronyResults2(Struc)

disp('finished')