
% Script to load the X matrix of beta values and perform the hyperalignment
% across subjects.


clearvars; close all;

baseDir = '~/Soundscapes';
path_results = fullfile(baseDir, 'results_balanced');
addpath(genpath(fullfile(baseDir, 'code/toolbox')));

hemi = {'lh' 'rh'};

% Define the ROIs => copy the gifti function
ROIsFolder = fullfile(baseDir, 'freesurfer_subject/fsaverage/atlasmgz/');

roiLabel = readtable(fullfile(ROIsFolder, 'Glasser2016_ColorLUT.txt'));
roiIdx = get_roi('fsaverage', ROIsFolder); % two columnn vector, one for left hemi and one for right hemi
list = roiLabel.LabelName_(2:end);

% list of subjects
subjects = {'sc001' 'sc002' 'sc003' 'sc004' 'sc005' 'sc006' 'sc007' 'sc008' 'sc009' 'sc010' ...
    'sc011' 'sc012' 'sc013' 'sc014' 'sc015' 'sc016' 'sc017' 'sc018' 'sc019' 'sc020'};

object_names = {'hum_Art','hum_Eet','hum_Lot','hum_Mar','hum_Poi','hum_Tyt', ...
    'ins_Art', 'ins_eet','ins_Lot','ins_Mar','ins_Poi','ins_Tyt', ...
    'ani_Api','ani_Hus','ani_Hyl','ani_Kis','ani_Tip','ani_Val'};

iROI = 1:180;
nExperiments = 5; 


%% Load the data aligned across subjects

for hemi_i = 1:length(hemi)
    
    for ROI_i = 1:length(iROI)
        
        % load the matrix of beta values (18 conditions * number of vertices)
        resultsDir = fullfile(path_results, strtrim(list{iROI(ROI_i)}), hemi{hemi_i}, 'auditory_objects');
        load(fullfile(resultsDir, 'data.mat'))
        
        disp([strtrim(list{iROI(ROI_i)}) ' ' hemi{hemi_i}])

        % Perform the hyperalignment across subjects, for each experiment
        aligned = cell(nExperiments, length(subjects)); 
        
        for exp_i = 1:nExperiments
            
            disp(['Exp ' num2str(exp_i)]);
            
            % Extract the data for the current experiment (1x5 cell)
            data_exp = data(exp_i, :);
            
            % Transpose each matrix in mean-centered data
            transposed_data = cellfun(@(x) x', data_exp, 'UniformOutput', false);  

            % Mean-center each row (subtract the mean across columns for each row)
            mean_centered_data = cellfun(@(x) x - mean(x, 2), transposed_data, 'UniformOutput', false);

            % Hyperalign the data across the subjects (assumes hyperalign accepts 1x5 cell input)
            [aligned_exp, transforms] = hyperalign_noScaling(mean_centered_data{:});  % Hyperalign across subjects  
            
            % Transpose the aligned data for each subject (700x18 to 18x700)
            for subj_i = 1:length(subjects)
                aligned{exp_i, subj_i} = aligned_exp{subj_i}';  % Transpose each subject's aligned data
            end 
            clear data_exp
    
        end
        
        save(fullfile(resultsDir, 'aligned_data.mat'), 'aligned')

    end
               
end



