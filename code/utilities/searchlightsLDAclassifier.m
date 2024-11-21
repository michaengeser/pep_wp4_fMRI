function searchlightsLDAclassifier(subs, nRuns, nTrials, map)
%
% The data used here is available from http://cosmomvpa.org/datadb.zip
%
% This example uses the following dataset:
% + 'digit'
%    A participant made finger pressed with the index and middle finger of
%    the right hand during 4 runs in an fMRI study. Each run was divided in
%    4 blocks with presses of each finger and analyzed with the GLM,
%    resulting in 2*4*4=32 t-values
%
% #   For CoSMoMVPA's copyright information and license terms,   #
% #   see the COPYING file distributed with CoSMoMVPA.           #


% reset citation list
cosmo_check_external('-tic');

%% LDA classifier searchlight analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This analysis identified brain regions where the categories can be
% distinguished using an odd-even partitioning scheme and a Linear
% Discriminant Analysis (LDA) classifier.

% loop through subjects
for iSub = 1:length(subs)

    if subs(iSub) < 10
        subID = ['sub-00', num2str(subs(iSub))];
    elseif subs(iSub) < 100
        subID = ['sub-0', num2str(subs(iSub))];
    end

    disp(['Start searchlight for subject ',  num2str(subs(iSub)), ' on ',...
        map, '-map']);

    if strcmp(map, 't')
        % Initialize dataset cell array
        datasets = cell(numel(nRuns), 1);

        % get mask
        % Load all masks
        mask_files = {};
        for mask_run = 1:nRuns
            mask_files{mask_run} = fullfile(pwd, '..', 'derivatives', subID, 'exp_glm1_norm', ...
                sprintf('run-%02d', mask_run), 'mask.nii');
        end
        intersection_mask = spm_read_vols(spm_vol(mask_files{1})) > 0;
        for i = 2:numel(mask_files)
            current_mask = spm_read_vols(spm_vol(mask_files{i})) > 0;
            intersection_mask = intersection_mask & current_mask;
        end
        % Save the intersection mask
        hdr = spm_vol(mask_files{1});
        intersection_mask_dir = fullfile(pwd, '..', 'derivatives', subID,...
            'exp_glm1_norm', 'intersection_mask.nii');
        hdr.fname = intersection_mask_dir;
        spm_write_vol(hdr, intersection_mask);


        % Loop over runs to load data
        counter = 0;
        for run = 1:nRuns

            % get dir for folder
            con_dir = fullfile(pwd, '..', 'derivatives', subID, 'exp_glm1_norm', ...
                sprintf('run-%02d',run));

            for trial = 1:nTrials
                counter = counter + 1;

                % Load the contrast image for this run and trial
                con_file = fullfile(con_dir, sprintf('con_%04d.nii',trial));

                % Convert to CoSMoMVPA dataset
                if trial <= nTrials/2; target = 1; else; target = 2;end
                ds = cosmo_fmri_dataset(con_file, ...
                    'mask', intersection_mask_dir, ... % Set brain mask
                    'targets', target, ... % Set condition labels (1 = bathroom)
                    'chunks', run);     % Set run identifiers

                % Store dataset
                datasets{counter} = ds;
            end
        end

        % Combine all runs into a single dataset
        ds_per_run = cosmo_stack(datasets);

    elseif strcmp(map, 'b')

        % get path
        betaPath = fullfile(pwd, '..', 'derivatives', subID, 'GLMsingle_betas', ...
            'beta.nii');

        % load beta map
        ds_per_run = cosmo_fmri_dataset(betaPath);

        % add chunks and targets
        nSamples = height(ds_per_run.samples);
        nSamplesPerRun = nSamples/nRuns;
        ds_per_run.sa.targets = repmat([ones(1,nTrials/2), ones(1,nTrials/2)*2,...
            ones(1,nSamplesPerRun-nTrials)*3], 1, nRuns)';
        preChunks = repmat(1:nRuns, nSamplesPerRun, 1);
        ds_per_run.sa.chunks = reshape(preChunks,[],1);

        % remove living room trials trials (target == 3)
        ds_per_run = cosmo_slice(ds_per_run, (ds_per_run.sa.targets ~= 3), 1); 

    else
        error('map not defined')
    end

    % Use the cosmo_cross_validation_measure and set its parameters
    % (classifier and partitions) in a measure_args struct.
    measure = @cosmo_crossvalidation_measure;
    measure_args = struct();

    % Define which classifier to use, using a function handle.
    % Alternatives are @cosmo_classify_{svm,matlabsvm,libsvm,nn,naive_bayes}
    measure_args.classifier = @cosmo_classify_lda;

    % Set partition scheme. odd_even is fast; for publication-quality analysis
    % nfold_partitioner is recommended.
    % Alternatives are:
    % - cosmo_nfold_partitioner    (take-one-chunk-out crossvalidation)
    % - cosmo_nchoosek_partitioner (take-K-chunks-out  "             ").
    measure_args.partitions = cosmo_oddeven_partitioner(ds_per_run);

    % print measure and arguments
    fprintf('Searchlight measure:\n');
    cosmo_disp(measure);
    fprintf('Searchlight measure arguments:\n');
    cosmo_disp(measure_args);

    % Define a neighborhood with approximately 100 voxels in each searchlight.
    nvoxels_per_searchlight = 500;
    nbrhood=cosmo_spherical_neighborhood(ds_per_run,...
        'count',nvoxels_per_searchlight);

    % Run the searchlight
    lda_results = cosmo_searchlight(ds_per_run,nbrhood,measure,measure_args);

    % Define output location
    output_fn=fullfile(con_dir, '..', '..', 'lda_searchlight', filesep);
    if ~exist(output_fn, 'dir')
        mkdir(output_fn);
    end

    % Store results to disc
    cosmo_map2fmri(lda_results, fullfile(output_fn, 'lda_searchlight.nii'));

end

%% average across participants

output_dir = fullfile(pwd, '..', 'derivatives', 'averaged_results');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Loop through each subject
for iSub = 1:length(subs)

    if subs(iSub) < 10
        subID = ['sub-00', num2str(subs(iSub))];
    elseif subs(iSub) < 100
        subID = ['sub-0', num2str(subs(iSub))];
    end

    % get subjects data
    decoding_map_file = fullfile(pwd, '..', 'derivatives', subID,...
        'lda_searchlight', 'lda_searchlight.nii');

    % Load decoding results using CoSMoMVPA
    if exist(decoding_map_file, 'file')
        ds = cosmo_fmri_dataset(decoding_map_file);

        % Store results for averaging
        all_results(iSub, :) = ds.samples; % Initialize with first subject's data
    else
        warning('Decoding map not found for %s, skipping...', subject_dirs(i).name);
    end
end

% take mean
average_results = mean(all_results);

% Save the averaged decoding map
fprintf('Saving averaged decoding map...\n');
template_ds = cosmo_fmri_dataset(decoding_map_file);
template_ds.samples = average_results;
average_output_file = fullfile(output_dir, 'average_decoding_map.nii');
cosmo_map2fmri(template_ds, average_output_file);

fprintf('Averaged decoding map saved to: %s\n', average_output_file);
%% second level 

% Set paths
second_level_dir = fullfile(pwd, '..', 'derivatives', 'second_level_analysis');
if ~exist(second_level_dir, 'dir')
    mkdir(second_level_dir);
end

% Define chance level (adjust for your classification task)
chance_level = 0.5; % Binary classification

% Collect subject maps
subject_maps = cell(length(subs), 1);

for i = 1:length(subs)

    if subs(iSub) < 10
        subID = ['sub-00', num2str(subs(iSub))];
    elseif subs(iSub) < 100
        subID = ['sub-0', num2str(subs(iSub))];
    end

    % get subjects data
    decoding_map_file = fullfile(pwd, '..', 'derivatives', subID,...
        'lda_searchlight', 'lda_searchlight.nii');
    
    if exist(decoding_map_file, 'file')
        subject_maps{i} = decoding_map_file;
    else
        warning('Missing decoding map for %s, skipping...', subject_dirs(i).name);
    end
end


% Create SPM batch for second-level analysis
spm('defaults', 'FMRI');
spm_jobman('initcfg');

matlabbatch = {};

% 1. Specify second-level model
matlabbatch{1}.spm.stats.factorial_design.dir = {second_level_dir};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = subject_maps;
matlabbatch{1}.spm.stats.factorial_design.cov = struct([]);
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct([]);
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1; % Implicit mask
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

% 2. Estimate the model
matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(second_level_dir, 'SPM.mat')};

% 3. Perform one-sample t-test
matlabbatch{3}.spm.stats.con.spmmat = {fullfile(second_level_dir, 'SPM.mat')};
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Above chance';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

% Run batch
spm_jobman('run', matlabbatch);

disp('Second-level analysis complete!');



end