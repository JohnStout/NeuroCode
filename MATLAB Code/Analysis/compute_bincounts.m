%% Compute_bincounts
% this function uses the freedman_diaconis rule and averages across animals
% to determine an average bin count 
%
% Written by John Stout
% last edit 12/11/18

%% fun
function [nbins_av] = compute_bincounts(input)

% define which part of the maze
    mazei = [1 6]; % 1 is stem entry 6 is choice-point exit

% Initialize variables to flip over all folders    
    if input.Prelimbic == 1;
        Datafolders = 'X:\01.Experiments\John n Andrew\Dual optogenetics w mPFC recordings\All Subjects - DNMP\Good performance\Prelimbic';
    elseif input.OFC ==1;
        Datafolders = 'X:\01.Experiments\John n Andrew\Dual optogenetics w mPFC recordings\All Subjects - DNMP\Good performance\Orbital Frontal';    
    elseif input.AnteriorCingulate == 1;
        Datafolders = 'X:\01.Experiments\John n Andrew\Dual optogenetics w mPFC recordings\All Subjects - DNMP\Good performance\Anterior Cingulate';
    elseif input.mPFC_good == 1;
        Datafolders = 'X:\01.Experiments\John n Andrew\Dual optogenetics w mPFC recordings\All Subjects - DNMP\Good performance\Medial Prefrontal Cortex';
    elseif input.mPFC_poor == 1;
        Datafolders = 'X:\01.Experiments\John n Andrew\Dual optogenetics w mPFC recordings\All Subjects - DNMP\Poor Performance'; 
    elseif input.VentralOrbital == 1;
        Datafolders = 'X:\01.Experiments\John n Andrew\Dual optogenetics w mPFC recordings\All Subjects - DNMP\Good performance\Ventral Orbital';
    elseif input.MedialOrbital == 1;
        Datafolders = 'X:\01.Experiments\John n Andrew\Dual optogenetics w mPFC recordings\All Subjects - DNMP\Good performance\Medial Orbital';
    else
        disp('Warning - Error in loading Datafolders')
    end
    
    % some last minute stuff
      cd(Datafolders);  
      folder_names = dir;
    
% flip across sessions    
for nn = 3:length(folder_names)
      Datafolders = Datafolders;
 
    % make a variable that keeps track of which session
      cd(Datafolders);
      folder_names = dir;
      temp_folder = folder_names(nn).name;
      cd(temp_folder);
      datafolder = pwd;
      cd(datafolder);
    
    % load int
    load(strcat(datafolder,'\Int_file.mat'));
    
    % load vt data
    load(strcat(datafolder, '\VT1.mat'));
    try load(strcat(datafolder,'\ExtractedYSample.mat'));
        %disp('Successfully loaded ExtractedY Sample')
    catch 
         %disp('ExtractedY from VT1 will be included');
    end
    
    % Do you want perfect sessions included?
    if input.notperfectsess == 1 
        B = any(Int(:,4));
        if B == 0
           Int = [];
                continue
        end
    end
    
% switch directories
cd 'X:\03. Lab Procedures and Protocols\MATLABToolbox\John code and edits\Information Theory'

% run freedman diaconis function
nbins = freedman_diaconis(Int,ExtractedY,TimeStamps_VT,mazei);

% store data
session_bins{nn} = nbins;

% house keeping
nbins = [];
    
end

% average across bins
nbins_av = ceil(mean(cell2mat(session_bins(3:end))));
