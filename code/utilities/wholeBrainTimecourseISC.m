function wholeBrain = wholeBrainTimecourseISC(d, cfg)


% evaluate input
if ~isfield(cfg, 'plotting'); cfg.plotting = true; end
if ~isfield(cfg, 'saving'); cfg.saving = false; end
if ~isfield(cfg, 'nRuns'); cfg.nRuns = 10; end
if ~isfield(cfg, 'regressOutMean'); cfg.regressOutMean = true; end
if ~isfield(cfg, 'brainMask'); cfg.brainMask = 'group_mask_thresholded'; end
if ~isfield(cfg, 'predictor_RDMs'); cfg.predictor_RDMs = {'typical_late', 'control_late'}; end
if ~isfield(cfg, 'RDM_to_partial_out'); cfg.RDM_to_partial_out = cfg.predictor_RDMs; end
if ~isfield(cfg, 'correlation_type'); cfg.correlation_type = 'pearson';end
if ~isfield(cfg, 'permutation_test'); cfg.permutation_test = false;end
if ~isfield(cfg, 'n_permutations'); cfg.n_permutations = 10000;end
if ~isfield(cfg, 'permutation_type'); cfg.permutation_type = 'row_col_shuffle_ref';end
if ~isfield(cfg, 'dnns'); cfg.dnns = {cfg.dnn};end

cfg.predictor_RDMs = {'typical_late', 'control_late'};
cfg.RDM_to_partial_out = cfg.predictor_RDMs;

% convert to local configurations
if cfg.plotting; plotting = true; else plotting = false; end
if cfg.saving; saving = true; else saving = false; end
cfg.plotting = false;
cfg.saving = false;

% get brain mask
makeBrainMask(cfg)
brainMaskPath = fullfile(pwd, '..', 'derivatives', 'group_level', [cfg.brainMask, '.nii']);
brainMask = load_untouch_nii(brainMaskPath);


% load brain mask
% loop through categories
for category = cfg.categories
    category = char(category);
    Category = strcat(upper(category(1)),lower(category(2:end))); % capitalize first letter

    % ISC path
    ISCpath = fullfile(pwd, '..', 'ISCtoolbox', category, 'memMaps.mat');
    load(ISCpath)

    % init 5D matrix (x, y, z, c = vectorized ISC mat, r = run)
    xyzcr = nan([size(memMaps.cormatMap.whole.band0.Session1.cor.Data.xyzc), cfg.nRuns]);
    xyzcr = single(xyzcr);

    % init predictor RDMs
    RDMs = [];
    RDMs = d.DNN.(cfg.dnn).control.(category).subject_mean(1); % placeholder RDM
    labels = {RDMs.name};
    [RDMs, cfg.labels] = evaluate_predictor_RDMs(d, RDMs, labels, cfg, category);
    cfg.labels{1} = 'VoxelTimecourse';
    for field = 1:numel({RDMs.name})
        RDMs(field).name = char(cfg.labels{field}); % give it a comprehensive name
    end

    % get mean across runs
    for iRun = 1:cfg.nRuns

        % get current matrix
        xyzcr(:,:,:,:,iRun) = single(memMaps.cormatMap.whole.band0.(['Session', num2str(iRun)]).cor.Data.xyzc);
    end
    meanXyzc = mean(xyzcr, 5, 'omitnan');

    % init correlation matrix
    reshapedMeanXyzc = reshape(meanXyzc, [], size(meanXyzc, 4));
    rVec = zeros(1, size(reshapedMeanXyzc, 1));
    pVec = ones(1, size(reshapedMeanXyzc, 1));

    parpool(8);
    RDMs(2).RDM(eye(cfg.n, cfg.n) == 1) =0;
    vec_RDM2 = squareform(RDMs(2).RDM);
    RDMs(3).RDM(eye(cfg.n, cfg.n) == 1) =0;
    vec_RDM3 = squareform(RDMs(3).RDM);
    parfor iVoxel = 1:size(reshapedMeanXyzc, 1)

        % skip if there is mean ISC = 0 for this voxel
        if mean(reshapedMeanXyzc(iVoxel, :)) == 0
            continue
        else

            % get RDM for this voxel
            vec_RDM1 = reshapedMeanXyzc(iVoxel, :);

            % partial correlation
            [rVec(iVoxel), pVec(iVoxel)] = partialcorr(vec_RDM1', vec_RDM2', vec_RDM3',...
                'Tail', 'right', 'Type', 'Pearson');

        end

        if mod(iVoxel/size(reshapedMeanXyzc, 1), 0.1) < 0.000001
            disp(['Progress: ' num2str(round((iVoxel/size(reshapedMeanXyzc, 1))*100)), '%'])
        end
        
    end

    % make 3D brain map again
    wholeBrain.(category).rMap = reshape(rVec, size(brainMask.img));
    wholeBrain.(category).SignificanceMap = reshape(pVec, size(brainMask.img));

    % write nifi file
    nii = brainMask;
    nii.hdr.dime.datatype = 16;
    nii.hdr.dime.bitpix = 32;
    nii.img = single(wholeBrain.(category).rMap);
    save_untouch_nii(nii, fullfile(pwd, '..', 'ISCtoolbox', category,...
        'wholeBrainCorMat.nii'));
    nii.img = single(wholeBrain.(category).SignificanceMap);
    save_untouch_nii(nii, fullfile(pwd, '..', 'ISCtoolbox', category,...
        'wholeBrainPvalMat.nii'));

    delete(gcp('nocreate'))
end

end
