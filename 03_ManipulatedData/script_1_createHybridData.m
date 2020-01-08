%% Script: Generate in silico manipulated data ("added spikes" and "deleted spikes") from source data
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

sourceFiles = dir([path filesep 'Data_Source']);

%% Load File

N=40; % number of trials per shuffled point
rec_dur=300;

%% Select Shuffling Methods:
shuffleMethods={ ...
    'addedNoise', ...
    'deletedSpikes_perElectrode', ...
    };
shuffleMethod=shuffleMethods(1:end);


%% for every Test-File
for f=3:size(sourceFiles,1)
    
    % load original data (=whole chip)
    disp(sourceFiles(f).name)
    temp=load([sourceFiles(f).folder filesep sourceFiles(f).name]);
    TS = temp.TS;
    TS(TS==0)=NaN;
    TS=sort(TS);
    clear temp
    
    x=0:0.1:1; % levels of manipulation (0% to 100%)
    
    %% for every shuffling method
    for m=1:size(shuffleMethod,2)
        
        %% for every percentage value
        for p=1:length(x)
            percentage=x(p);
            disp(num2str(percentage*100))
            %% for every sample
            for n=1:N
                clear M_TS
                
                % add random spikes (noise):
                if strcmp(shuffleMethod{m},'addedNoise')
                    TSrand=TS;
                    a=1/10000; b=rec_dur;
                    
                    for el=1:size(TS,2)
                        tmp = TS(:,el);
                        numSpikes=length(tmp(~isnan(tmp)));
                        numberOfSpikes = int32(percentage*numSpikes*1.0); % *0.1: 100% manipulation means, that 0.1*length(TS) of noise spikes has been added, so maximal 10 % of spikes are added (100% would make not a big difference anymore, sync measures are more sensitive from 0 to 10 % )
                        noise= a + (b-a).*rand(numberOfSpikes,1);
                        if ~isempty(noise) && numSpikes>5*rec_dur/60 % more than five spikes per minute (here five minutes) to consider electrode as active
                            for ii=1:length(noise)
                                TSrand(end+1,el)=noise(ii); % add noise to end of array
                                
                                while length(unique(TSrand(:,el))) ~= length(TSrand(:,el))
                                    noise2= a + (b-a).*rand(1,1);
                                    TSrand(end,el)=noise2;
                                    disp('same TS')
                                    len_TSrand=length(TSrand(:,el));
                                    len_unique=length(unique(TSrand(:,el)));
                                end
                            end
                            TSrand(TSrand==0)=NaN;
                            TSrand=sort(TSrand);
                        end
                    end
                    M_TS=TSrand; % rename matrix
                end
                
                
                
                % delete random spikes:
                if strcmp(shuffleMethod{m},'deletedSpikes_perElectrode')
                    
                    TSrand=TS;
                    TSrand(TSrand==0)=NaN;
                    TSrand=sort(TSrand);
                                    
                    Xmax = 0.9; % if percentage = 1 it should be 90 % (so 10 % of spikes remain)
                    Xmin = 0; % if percentage = 0 it should be 0 %
                    X = Xmin + ((percentage)*(Xmax-Xmin))/(1);     
                    
                    for el=1:size(TS,2)
                        tmp=TS(:,el);
                        numSpikes=length(tmp(~isnan(tmp)));
                        numberOfShuffles = int32(X*numSpikes);
                        numSpikesAfterDeletion = numSpikes - numberOfShuffles;
                        
                        if numSpikesAfterDeletion >= 6*rec_dur/60 % only delete spikes if at least 6 spikes per minute are on current electrode after spikes have been deleted
                             idxDelete = randperm(numSpikes,numberOfShuffles); % pick n=numberOfShuffles random numbers in the range from 1:numSpikes (instead of using randi, here randperm is used to avoid repeated integer numbers)
                             TSrand(idxDelete,el)=NaN; % delete spike
                        end
                    end
                    
                    TSrand=sort(TSrand);
                    M_TS=TSrand; % rename matrix
                end

                
                %% save files
                folder_name = ['Data_Shuffled' filesep num2str(f) '-' sourceFiles(f).name filesep num2str(m) '-' shuffleMethod{m} filesep num2str(p) '-percentage' filesep]; % ShuffledData/File#/ShuffleMethod/percentage
                if ~exist(folder_name,'dir')
                    mkdir(folder_name);
                end
                filename=num2str(n);
                S.M_TS = M_TS;
                S.percentage = percentage;
                S.percentage_index = p;
                save([folder_name filename], 'S')
                disp(['saved: ' folder_name filename])
                
            end % END: for all samples n
        end % END: for all percentages x
    end % END: for all shuffling methods
end % END: for all Test files
disp('finished')