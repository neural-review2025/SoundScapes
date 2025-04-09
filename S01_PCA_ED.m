
% Script to perform the PCA at subject-level, align the scoes with the procrustes alignment and 
% extract Principal angle and euclidian distances, fro each ROI;
% Perform statistical analysis to test for regions in which PA or ED
% significantly changes according to the soundscape.

clearvars; close all;

baseDir = '/Users/cp3488/Documents/Projects/Helsinki/Soundscapes/Conference/github/Soundscapes';
path_results = fullfile(baseDir, 'results_balanced');
addpath(genpath(fullfile(baseDir, 'code/toolbox')));
addpath(genpath(getenv('FREESURFER_HOME'))); % Ensure FREESURFER_HOME is set
setenv('SUBJECTS_DIR', fullfile(baseDir, 'freesurfer_subject'));

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

nDim = 3; % Number of dimensions to keep in the pca
nExperiments = 5; 
fontsize = 18;
nPermutations = 50000;  % Number of permutations for subject-level analysis

% Define category indices
nCategories = 3; % Humans, Instruments, Animals
categoryIdx = {
    1:6;   % Humans
    7:12;  % Instruments
    13:18  % Animals
    };

iROI_selected = [124 125 107 129 75 54];
%iROILabels = {'PBelt', 'A5', 'TA2', 'STSdp', '45', '6d'}';

label_exp = 'aligned_NoScaling';



%% Load the aligned data, run the PCA, extract Principal Angle and Euclidia Distance and test for significant differences at ROI-level

for hemi_i = 1:length(hemi)
    
    for ROI_i = 1:length(iROI)
        
        % load the X matrix after hyperalignment
        resultsDir = fullfile(path_results, strtrim(list{iROI(ROI_i)}), hemi{hemi_i}, 'auditory_objects');
        load(fullfile(resultsDir, 'aligned_data.mat'))
        
        disp([strtrim(list{iROI(ROI_i)}) ' ' hemi{hemi_i}])

        for sub_i = 1:length(subjects)
            
            % Perform the PCA
            X_exp = cell(1, nExperiments);
            for expIdx = 1:nExperiments
                X_exp{expIdx} = aligned{expIdx,sub_i};
                X_exp{expIdx} = (X_exp{expIdx} - mean(X_exp{expIdx}, 1));
                [~, scores{expIdx},~,~,explained{sub_i}(:,expIdx)] = pca(X_exp{expIdx});
            end
            
            % Align the scores across experiments
            for exp_i = 1:size(aligned,1)                
                numeric_data = scores;
                reference_score = numeric_data{3};  % Choose experiment 3 as reference
                if exp_i == 3
                    aligned_tmp{exp_i} = reference_score;
                    procrustes_distances(sub_i, exp_i) = 0;  % Reference experiment distance is 0
                else
                    [d, aligned_matrix_score, transform] = procrustes(reference_score, numeric_data{exp_i}, 'Scaling', false);
                    aligned_tmp{exp_i} = aligned_matrix_score;
                    procrustes_distances(sub_i, exp_i) = d; % Store Procrustes distance
                end                
            end
            Y_exp = aligned_tmp;
            
            % Select the experiments we want to analyze
            exp = [2,3,5];
            categoryNames = {'Humans', 'Instruments', 'Animals'};

            % Extract the centroid of each category
            for exp_ii = 1:length(exp)
                exp_i = exp(exp_ii);
                for catIdx = 1:nCategories
                    centroids_all{sub_i, exp_ii}(catIdx, :) = mean(Y_exp{exp_i}(categoryIdx{catIdx}, 1:nDim), 1); 
                end
            end
        
            % Compute Euclidean distances between category centroids
            for exp_ii = 1:length(exp)
                exp_i = exp(exp_ii);
                
                for catIdx = 1:nCategories           

                    % Compute all-to-all Euclidean distances
                    X2 = Y_exp{3}(categoryIdx{catIdx}, 1:nDim); % Always Experiment 3
                    X1 = Y_exp{exp_i}(categoryIdx{catIdx}, 1:nDim); % Current experiment (2 attended or 2 distarctors)

                    % Compute Euclidean distance for each corresponding row
                    distances = sqrt(sum((X1 - X2).^2, 2));
                    % Compute the mean Euclidean distance
                    mean_distance{sub_i}(catIdx, exp_ii) = mean(distances);
                    
                    % Calculate the absolute distance between each stimulus
                    % and the centroid, to plot the elipsoid around the centroid in figure 2.
                    stimuli = Y_exp{exp_i}(categoryIdx{catIdx}, 1:nDim);
                    for stimIdx = 1:size(stimuli, 1)
                        for pc_i = 1:nDim
                            % 6 rows * 3 PCs
                            abs_distances_within{exp_ii,catIdx}(stimIdx,pc_i) = sum(abs(stimuli(stimIdx, pc_i) - centroids_all{exp_ii}(catIdx, pc_i)));
                        end
                    end
                    
                    abs_distances_within_mean{sub_i, exp_ii}(catIdx, :)= mean(abs_distances_within{exp_ii,catIdx},1);
 
                end
            end
            
            % Run pca on the different categories in 3D in preparation for plane fitting
            eigvecs = cell(nCategories, nExperiments);
            % Compute PCA separately for each category and experiment
            for catIdx = 1:nCategories
                for expIdx = 1:nExperiments
                    [eigvecs{catIdx, expIdx},~,~,~, explainedVariance{catIdx, expIdx}] = pca(Y_exp{expIdx}(categoryIdx{catIdx}, 1:nDim));
                end
            end

            % Compute the Principal Angle
            for catIdx = 1:nCategories
                for exp_ii = 1:length(exp)
                    exp_i = exp(exp_ii);
                    cosTheta(catIdx, exp_ii) = planeAngle( ...
                        eigvecs{catIdx, exp_i}(:,1), eigvecs{catIdx, exp_i}(:,2), ...
                        eigvecs{catIdx, 3}(:,1), eigvecs{catIdx, 3}(:,2));     
                    theta{sub_i}(catIdx, exp_ii) = acosd(cosTheta(catIdx, exp_ii));
                end
            end 
            
        end
        
        
        %% Test for significant difference with permutations
        pValues{ROI_i, hemi_i} = zeros(nCategories, 2);  
        t_perm{ROI_i, hemi_i} = zeros(nCategories, 2);  
        t_value{ROI_i, hemi_i} = zeros(nCategories, 2);  

        for var_i = 1:2
            if var_i == 1
                x = theta; % Principal Angle
            else
                x = mean_distance; % ED using the pairwise distances
            end

            for catIdx = 1:nCategories
                for sub_i = 1:length(subjects)
                    exp1Att(sub_i,catIdx) = x{sub_i}(catIdx, 1);  
                    exp2Att(sub_i,catIdx) = x{sub_i}(catIdx, 3);
                end
                [h, p{ROI_i, hemi_i}, ci, stats] = ttest(exp1Att(:,catIdx), exp2Att(:,catIdx));
                t_value{ROI_i, hemi_i}(catIdx, var_i) = stats.tstat; 
                
                % Compute observed paired differences
                diffs = exp1Att(:,catIdx) - exp2Att(:,catIdx);
                observed_diff = mean(diffs);
                % Initialize a matrix to store permutation differences
                perm_diffs = zeros(nPermutations, 1);
                % Perform sign-flipping permutation test
                for permIdx = 1:nPermutations
                    % Randomly flip signs (+1 or -1)
                    flipSigns = randi([0,1], size(diffs)) * 2 - 1; 
                    permuted_diffs = diffs .* flipSigns;
                    % Compute permuted difference mean
                    perm_diffs(permIdx) = mean(permuted_diffs);
                end         
                % Compute permutation-based standard deviation
                perm_std = std(perm_diffs);

                % Compute permutation t-value
                t_perm{ROI_i, hemi_i}(catIdx, var_i) = observed_diff / perm_std;

                % Compute p-value (two-tailed test)
                pValues{ROI_i, hemi_i}(catIdx, var_i) = mean(abs(perm_diffs) >= abs(observed_diff));
            end

        end
        
        if any(ROI_i == iROI_selected)
                    
            figureDir = fullfile(path_results, 'results', label_exp, strtrim(list{iROI(ROI_i)}));
            mkdir(figureDir)
            
            % 1. Plot the centroids for each category and experiment
            figure('units', 'inches', 'position', [1 1 10 8]);
            exp_labels = {'3OAA', 'OA', '3OAD'};
            colors = {'r', 'g', 'b'};
            markers = {'o', 's', 'v'};

            % Create empty cell arrays to hold handles for the legend
            h_exp = cell(1, nExperiments);
            h_cat = cell(1, nCategories);

            for exp_i = 1:length(exp)
                hold on;
                centroid_mean = zeros(nCategories, 3);
                centroid_std = zeros(nCategories, 3);

                % Compute mean and std centroids across subjects
                for i = 1:nCategories
                    centroid_data = cell2mat(cellfun(@(x) x(i, :), centroids_all(:, exp_i), 'UniformOutput', false));
                    distance = cell2mat(cellfun(@(x) x(i, :), abs_distances_within_mean(:, exp_i), 'UniformOutput', false));
                    centroid_mean(i, :) = mean(centroid_data, 1);
                    centroid_sem(i, :) = std(distance, 0, 1) / sqrt(size(distance, 1));
                end

                for i = 1:nCategories
                    % Plot mean centroids and store handle for legend
                    h_cat{i} = plot3(centroid_mean(i, 1), centroid_mean(i, 2), centroid_mean(i, 3), ...
                        markers{i}, 'MarkerSize', 10, 'MarkerFaceColor', colors{exp_i}, ...
                        'DisplayName', categoryNames{i});
                    h_cat{i} = plot3(nan, nan, nan, markers{i}, 'Color', 'k', 'MarkerFaceColor', 'k','DisplayName', categoryNames{i});

                end
                h_exp{exp_i} = plot3(nan, nan, nan, colors{exp_i}, 'LineWidth', 2, 'DisplayName', exp_labels{exp_i});

                % Draw ellipsoids representing SEM for each category
                for i = 1:nCategories
                    [x, y, z] = ellipsoid(centroid_mean(i, 1), centroid_mean(i, 2), centroid_mean(i, 3), ...
                                           centroid_sem(i, 1), centroid_sem(i, 2), centroid_sem(i, 3), 10);
                    surf(x, y, z, 'FaceColor', colors{exp_i}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                end

                % Draw lines between mean centroids
                for i = 1:nCategories
                    for j = i+1:nCategories
                        mid_x = (centroid_mean(i,1) + centroid_mean(j,1)) / 2;
                        mid_y = (centroid_mean(i,2) + centroid_mean(j,2)) / 2;
                        mid_z = (centroid_mean(i,3) + centroid_mean(j,3)) / 2;

                        plot3([centroid_mean(i,1), centroid_mean(j,1)], ...
                            [centroid_mean(i,2), centroid_mean(j,2)], ...
                            [centroid_mean(i,3), centroid_mean(j,3)], [colors{exp_i} '-'], 'LineWidth', 1.5);

                    end
                end

                ax = gca; % Get current axes
                ax.FontSize = 16; % Adjust as needed
                grid off; % Remove internal grid lines
                box on; % Make the outer cube visible
                axis square;

                % Increase the thickness of the box lines
                ax.XAxis.LineWidth = 2; % Adjust the thickness for the X axis line
                ax.YAxis.LineWidth = 2; % Adjust the thickness for the Y axis line
                ax.ZAxis.LineWidth = 2; % Adjust the thickness for the Z axis line

                xlabel('PC1', 'FontSize', fontsize);
                ylabel('PC2', 'FontSize', fontsize);
                zlabel('PC3', 'FontSize', fontsize);
                rotate3d on;
                set(gca, 'XTickLabel', [], 'YTickLabel', [], 'ZTickLabel', []);
            end

            % Optional: save the figures
            saveas(gcf, fullfile(figureDir, ['ED_Elipsoids_' strtrim(list{iROI(ROI_i)}) '_' hemi{hemi_i} '.jpg']));
            exportgraphics(gcf, fullfile(figureDir, ['ED_Elipsoids_' strtrim(list{iROI(ROI_i)}) '_' hemi{hemi_i} '.pdf']), 'ContentType', 'vector');

            % 2. Plot the values of ED
            ED_data = [];  % To store the theta values
            categoryNames = {'hum', 'ins', 'ani'};

            % Loop through categories and experiments
            for catIdx = 1:nCategories
                for expIdx = 1:length(exp) % Using only the first two experiments for comparison  
                     if expIdx == 2
                     else
                        % Loop through subjects and extract the corresponding theta values
                        for sub_i = 1:length(subjects)
                            ED_values(sub_i) = squeeze(mean_distance{sub_i}(catIdx, expIdx)); % ED based on the pairwise distances
                        end

                        % Reshape theta_values to a column vector and add to the data array
                        ED_data = [ED_data ED_values'];
                    end
                end
            end

            data_to_plot = ED_data(:, 1:2); % choose the category -> 1:2 humand; 3:4 instruments; 5:6 animals
            % Calculate the mean for each condition (across subjects)
            means = mean(data_to_plot, 1);
            % Calculate the standard deviation for each condition (across subjects)
            stds = std(data_to_plot, 0, 1);
            % Calculate the number of subjects (rows in data_to_plot)
            N = size(data_to_plot, 1);
            % Calculate the Standard Error of the Mean (SEM)
            sem = stds / sqrt(N);  % SEM = std / sqrt(N)

            % Create the plot
            figure;
            hold on;
            % Plot individual subject trends (lines for each subject in grey)
            for i = 1:size(data_to_plot, 1)  % Loop over 20 subjects
                plot(1:2, data_to_plot(i,:), 'Color', [0.5 0.5 0.5], 'LineWidth', 4);  % Plot each subject's data in grey
            end
            % Plot the mean for each condition (blue circles)
            plot(1:2, means, 'ro-', 'MarkerFaceColor', 'r', 'MarkerSize', 1, 'LineWidth', 10);
            % Plot the SEM as error bars
            errorbar(1:2, means, sem, 'r', 'LineWidth', 10, 'CapSize', 10);
            % Set x-axis limits to go from 0 to 3
            xlim([0 3]);
            % Customize the plot
            xlabel('Condition', 'FontSize', 14);
            ylabel('Euclidean Distance', 'FontSize', 14);
            xticks(1:2);
            xticklabels({'Condition 1', 'Condition 2'});  % Customize based on your labels
            xtickangle(45);  % Rotate x-axis labels for readability
            hold off;
            
        end
        
        close all;
 
    end
              
end

%% Save the statistics
StatDir = fullfile(path_results, 'results', label_exp);
save(fullfile(StatDir, 'pValues.mat'), 'pValues') 
save(fullfile(StatDir, 't_perm.mat'), 't_perm') 
save(fullfile(StatDir, 't_value.mat'), 't_value') 

% Extract significant results after FDR-correction
categoryNames = {'Hum', 'Ins', 'Ani'};
repeatedCategories = repmat(categoryNames, 1, max(iROI))';
repeatedLabels = repelem(list, size(categoryNames,2));

load(fullfile(StatDir, 'pValues.mat'))
% Flatten all p-values into a single vector
p_all = cell2mat(pValues); % Convert cell array to numeric matrix
results = [repeatedLabels, repeatedCategories, num2cell(p_all)];
save(fullfile(StatDir, 'pValues_table.mat'), 'results') 

p_all = p_all(:); % Flatten into a column vector
% Perform FDR correction using mafdr
q_all = mafdr(p_all, 'BHFDR', true); 

% Reshape corrected p-values back to the original 3x2 structure in each cell
q_matrix = reshape(q_all, [], 4);
q_cells = num2cell(q_matrix);
results_fdr = [repeatedLabels, repeatedCategories, q_cells];
save(fullfile(StatDir, 'pValues_fdr.mat'), 'results_fdr') 



%% Plot the results on the brain surface
% choose the category to plot

% Loop over each category
for categoryIdx = 1:length(categoryNames)
    category = categoryNames{categoryIdx};
    
    % Load p-values and t-values for the current category
    load(fullfile(StatDir, 'pValues_fdr.mat'))
    p_values = results_fdr(strcmp(results_fdr(:,2), category), :);
    p_values_lh = cell2mat(p_values(:,4)); % Extract FDR-corrected p-values
    p_values_rh = cell2mat(p_values(:,6)); % Extract FDR-corrected p-values
    p_values_all = [p_values_lh p_values_rh];

    load(fullfile(StatDir, 't_perm.mat'))
    repeatedCategories = repmat(categoryNames, 1, max(iROI))';
    repeatedLabels = repelem(list, size(categoryNames,2));

    % Flatten all p-values into a single vector
    t_perm = cell2mat(t_perm); % Convert cell array to numeric matrix
    results = [repeatedLabels, repeatedCategories, num2cell(t_perm)];
    t_values = results(strcmp(results(:,2), category), :);
    t_values_lh = cell2mat(t_values(:,4)); % Extract FDR-corrected p-values
    t_values_rh = cell2mat(t_values(:,6)); % Extract FDR-corrected p-values
    t_values_all = [t_values_lh t_values_rh];

    hSize = get_surfsize('fsaverage'); 
    hSizeIdx = [1, hSize(1); hSize(1)+1, sum(hSize)]; % Hemisphere indexing

    % Initialize p-value maps with NaN (or zeros) to store values for all ROIs
    p_map_lh = NaN(size(roiIdx, 1), 1);  % Left hemisphere
    p_map_rh = NaN(size(roiIdx, 1), 1);  % Right hemisphere
    t_map_lh = NaN(size(roiIdx, 1), 1);  % Left hemisphere
    t_map_rh = NaN(size(roiIdx, 1), 1);  % Right hemisphere

    for ROI_i = 1:length(iROI)
        % Get the indices of vertices belonging to each ROI in the Glasser Atlas
        indx_lh = (roiIdx(:,1) == iROI(ROI_i));
        indx_rh = (roiIdx(:,2) == iROI(ROI_i));

        % Assign the p-value to all vertices in the ROI (keeping previous values)
        p_map_lh(indx_lh) = p_values_all(ROI_i,1);
        p_map_rh(indx_rh) = p_values_all(ROI_i,2);
        t_map_lh(indx_lh) = t_values_all(ROI_i,1);
        t_map_rh(indx_rh) = t_values_all(ROI_i,2);
    end

    % Plot the results on the lateral surface of the brain
    vals = [t_map_lh; t_map_rh]; 
    vals_p = [p_map_lh; p_map_rh]; 
    vals(vals_p>0.05) = nan;

    [~,~, rgbimg] = cvnlookup('fsaverage',6,vals,[-4 -1],flipud(jet),[],[],[],{'rgbnan',1, 'roiname',{'Glasser2016'},'roicolor',{'k'},'drawroinames', 0});

    % Plot and save the image
    a = imagesc(rgbimg); axis image tight; axis off;
    % Setup color bar
    colormap(flipud(jet)) % color map
    caxis([-4 -1]) % range of the colorbar    
    cb = colorbar('SouthOutside');
    set(gca,'FontSize',20)
    ylabel(cb,'t-values','Rotation',360);
    
    % Export the plot for each category
    exportgraphics(gca, fullfile(StatDir, [category '_tvalues.jpg']))
    exportgraphics(gca, fullfile(StatDir, [category '_tvalues.pdf']))
end
