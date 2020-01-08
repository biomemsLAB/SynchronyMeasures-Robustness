%% Script: plot spike trains (raster plot) of the experimental data (without BIC and with 10 µM BIC)  
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
path_pics = [path filesep 'MEA-TS-Data' filesep 'MC_Th_Factor_5' filesep 'Pics' filesep 'Rasterplots'];


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
       end
   end
end


%% plot
if 1      
    vm = [0.1 0.2 0.5 1 2 5 10]; % bic concentrations
    Rows = size(vm,2);    
    allChips = {'18812','19498','19579','20692','24370','24659','25050','25165'};
    Cols=size(allChips,2); % number of chips

    row = 1;
    col=1;
    
    for f_idx=1:10 % index of file that is plotted, 1...5: ref, 6...10: bic
    
        hf=figure('Name','1','visible','off');

        for e=1:size(con,2) % for every concentration

            for c=1:size(con(e).chip,2) % for every Chip

                [~,col]=max(strcmp(allChips,con(e).chip(c).name));
                row=e;

                hs = subplot(Rows, Cols, Cols*(row-1)+col);
                hs=plotSpikeTrain(con(e).chip(c).file(f_idx).SPIKEZ.TS, hs, 'dot', 1);
                hs.XLabel.String=[num2str(vm(e)) 'uM BIC'];
                hs.Title.String=con(e).chip(c).name;

            end

        end

        if 1
            if ~exist(path_pics,'dir')
                    mkdir(path_pics);  
            end
                   
            hf.PaperPositionMode = 'auto';
            hf.Units='points';
            hf.Position=[0,0,1920,1080];
            print([path_pics filesep 'TS_f_idx_' num2str(f_idx)],'-dpng','-r600')
        end
    
    end

end
disp('finished')
