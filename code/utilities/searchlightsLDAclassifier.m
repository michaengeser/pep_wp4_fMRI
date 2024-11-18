function searchlightsLDAclassifier(subs, nRuns, nTrials)
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

    % Sanity check the dataset
    if cosmo_check_dataset(ds_per_run)
        disp([subID, ' - dataset checked'])
    else
        warning([subID, ' - dataset failed check'])
    end

    % print dataset
    fprintf('Dataset input:\n');
    cosmo_disp(ds_per_run);

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
    nvoxels_per_searchlight=100;
    nbrhood=cosmo_spherical_neighborhood(ds_per_run,...
        'count',nvoxels_per_searchlight);

    % Run the searchlight
    lda_results = cosmo_searchlight(ds_per_run,nbrhood,measure,measure_args);

    % print output dataset
    fprintf('Dataset output:\n');
    cosmo_disp(lda_results);

    % Plot the output
    cosmo_plot_slices(lda_results);

    % Define output location
    output_fn=fullfile(con_dir, '..', '..', 'lda_searchlight', filesep);
    if ~exist(output_fn, 'dir')
        mkdir(output_fn);
    end

    % Store results to disc
    cosmo_map2fmri(lda_results, fullfile(output_fn, 'lda_searchlight.nii'));

    % Show citation information
    cosmo_check_external('-cite');
end
end