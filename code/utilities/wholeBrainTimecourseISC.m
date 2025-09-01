function d = wholeBrainTimecourseISC(d, cfg)
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
if ~isfield(cfg, 'partial_correlation_type'); cfg.partial_correlation_type = 'Pearson';end
if ~isfield(cfg, 'dnns'); cfg.dnns = {cfg.dnn};end
% generate permutated subjects list
rng('default'); rng(1) % ensure reproducible outcome
cfg.random_seqs = cell(1, cfg.n_permutations);
RDM_template = squareform(1:nchoosek(cfg.n, 2));
for i = 1:cfg.n_permutations
    cfg.random_seqs{i} = randperm(cfg.n);
    RDM_shuffle = RDM_template(cfg.random_seqs{i}, cfg.random_seqs{i});
    cfg.random_RDM_vecs{i} = squareform(RDM_shuffle);
end
partial_correlation_type = cfg.partial_correlation_type;

% make for both categories parallel
warning('load both categories and evaluate at the same time so average can be computed on the fly')


cfg.RDM_to_partial_out = cfg.predictor_RDMs;

% convert to local configurations
if cfg.plotting; plotting = true; else plotting = false; end
if cfg.saving; saving = true; else saving = false; end
cfg.plotting = false;
cfg.saving = false;

% get brain mask
brainMaskPath = fullfile(pwd, '..', 'derivatives', 'group_level', [cfg.brainMask, '.nii']);
if ~exist(brainMaskPath, 'file')
    makeBrainMask(cfg)
end
brainMask = load_untouch_nii(brainMaskPath);


% load brain mask
% loop through categories
for category = cfg.categories
    category = char(category);
    Category = strcat(upper(category(1)),lower(category(2:end))); % capitalize first letter

    % ISC path
    % ISC path
    ISCpath = fullfile(pwd, '..', 'ISCtoolbox', category, 'results', 'memMaps.mat');
    load(ISCpath)

    % init 5D matrix (x, y, z, c = vectorized ISC mat, r = run)
    xyzcr = nan([size(memMaps.cormatMap.whole.band0.Session1.cor.Data.xyzc), cfg.nRuns]);

    % get mean across runs
    for iRun = 1:cfg.nRuns/2

        % get current matrix
        xyzcr(:,:,:,:,iRun) = single(memMaps.cormatMap.whole.band0.(['Session', num2str(iRun)]).cor.Data.xyzc);
    end
    xyzc = mean(xyzcr, 5, 'omitnan');
    xyzcr = [];

    % init correlation matrix
    reshapedxyzc.(category) = reshape(xyzc, [], size(xyzc, 4));
    xyzc = [];

end

% prepare variables for parallel computing
isc_values_bat = reshapedxyzc.bathroom;
isc_values_kit = reshapedxyzc.kitchen;

%% init predictor RDMs and get residuals of predictors
% (so we can run simple correlations in the loop below)

% bathroom
RDMs = d.DNN.(cfg.dnn).control.bathroom.subject_mean(1); % placeholder RDM
labels = {RDMs.name};
[RDMs, cfg.labels] = evaluate_predictor_RDMs(d, RDMs, labels, cfg, 'bathroom');

for iPred = 1:length(cfg.predictor_RDMs)
    RDMs(iPred+1).RDM(eye(cfg.n, cfg.n) == 1) = 0;
    bat_preds(:, iPred) = squareform(RDMs(iPred+1).RDM);
end

y_resid = regress(bat_preds(:, 1), [bat_preds(:, 2:end), ones(size(bat_preds(:, 2:end),1),1)]);
y_hat = [bat_preds(:, 2:end), ones(size(bat_preds(:, 2:end),1),1)] * y_resid;
bat_resid = bat_preds(:, 1) - y_hat;

% kitchen
RDMs = d.DNN.(cfg.dnn).control.kitchen.subject_mean(1); % placeholder RDM
labels = {RDMs.name};
[RDMs, cfg.labels] = evaluate_predictor_RDMs(d, RDMs, labels, cfg, 'kitchen');

for iPred = 1:length(cfg.predictor_RDMs)
    RDMs(iPred+1).RDM(eye(cfg.n, cfg.n) == 1) = 0;
    kit_preds(:, iPred) = squareform(RDMs(iPred+1).RDM);
end

y_resid = regress(kit_preds(:, 1), [kit_preds(:, 2:end), ones(size(kit_preds(:, 2:end),1),1)]);
y_hat = [kit_preds(:, 2:end), ones(size(kit_preds(:, 2:end),1),1)] * y_resid;
kit_resid = kit_preds(:, 1) - y_hat;

% if isempty(gcp('nocreate'))
%     parpool(10);
% end
for iVoxel = 1:size(isc_values_bat, 1)

    % make nan if mean ISC = 0 for this voxel
    if mean(isc_values_bat(iVoxel, :)) == 0 || ...
            mean(isc_values_kit(iVoxel, :)) == 0
        isc_values_bat(iVoxel, :) = nan;
        isc_values_kit(iVoxel, :) = nan;
    end
end

%% get correlation maps and save them

% write nifi file
nii = brainMask;
nii.hdr.dime.datatype = 16;
nii.hdr.dime.bitpix = 32;

% bathroom
bat_map = corr(isc_values_bat', bat_resid, 'Tail', 'right', 'Type', partial_correlation_type);
nii.img = single(reshape(bat_map, size(brainMask.img)));
save_untouch_nii(nii, fullfile(pwd, '..', 'ISCtoolbox', 'bathroom',...
    'wholeBrainCorMap.nii'));

% kitchen
kit_map = corr(isc_values_kit', kit_resid, 'Tail', 'right', 'Type', partial_correlation_type);
nii.img = single(reshape(kit_map, size(brainMask.img)));
save_untouch_nii(nii, fullfile(pwd, '..', 'ISCtoolbox', 'kitchen',...
    'wholeBrainCorMap.nii'));

% average
avg_map = (bat_map + kit_map)/2;
nii.img = single(reshape(avg_map, size(brainMask.img)));
save_untouch_nii(nii, fullfile(pwd, '..', 'ISCtoolbox', ...
    'wholeBrainCorMapAverage.nii'));


%% permutation test
if cfg.permutation_test
    % preallocate
    nHeigher_bat = zeros(length(bat_map), 1);
    nHeigher_bat(isnan(bat_map)) = nan;
    nHeigher_kit = zeros(length(kit_map), 1);
    nHeigher_kit(isnan(kit_map)) = nan;
    maxR_bat = zeros(cfg.nPermutations, 1);
    maxR_kit = zeros(cfg.nPermutations, 1);
    nHeigher_avg = nHeigher_bat;
    maxR_avg = maxR_bat;
    random_RDM_vecs = cfg.random_RDM_vecs;

    for perm = 1:cfg.n_permutations

        % bathroom
        perm_bat_map = corr(isc_values_bat(:, random_RDM_vecs{perm})', bat_resid,...
            'Tail', 'right', 'Type', partial_correlation_type);

        % test if heigher than actual correlation value
        heigherR = (perm_bat_map >= bat_map);
        nHeigher_bat = nHeigher_bat + heigherR;

        % store max value
        maxR_bat(perm) = max(perm_bat_map);


        % kitchen
        perm_kit_map = corr(isc_values_kit(:, random_RDM_vecs{perm})', kit_resid,...
            'Tail', 'right', 'Type', partial_correlation_type);

        % test if heigher than actual correlation value
        heigherR = (perm_kit_map >= kit_map);
        nHeigher_kit = nHeigher_kit + heigherR;

        % store max value
        maxR_kit(perm) = max(perm_kit_map);

        % average
        perm_avg_map = (perm_bat_map + perm_kit_map)/2;

        % test if heigher than actual correlation value
        heigherR = (perm_avg_map >= avg_map);
        nHeigher_avg = nHeigher_avg + heigherR;

        % store max value
        maxR_avg(perm) = max(perm_avg_map);
    end

    % get p values
    pMap_bat = nHeigher_bat/cfg.nPermutations;
    nii.img = single(reshape(pMap_bat, size(brainMask.img)));
    save_untouch_nii(nii, fullfile(pwd, '..', 'ISCtoolbox', 'bathroom',...
        'wholeBrainPvals.nii'));

    pMap_kit = nHeigher_kit/cfg.nPermutations;
    nii.img = single(reshape(pMap_kit, size(brainMask.img)));
    save_untouch_nii(nii, fullfile(pwd, '..', 'ISCtoolbox', 'bathroom',...
        'wholeBrainPvals.nii'));

    pMap_avg = nHeigher_avg/cfg.nPermutations;
    nii.img = single(reshape(pMap_avg, size(brainMask.img)));
    save_untouch_nii(nii, fullfile(pwd, '..', 'ISCtoolbox', ...
        'wholeBrainPvalsAverage.nii'));

    % get p values FDR-corrected
    pMap_bat(isnan(bat_map)) = nan;
    [~, ~, ~, pMap_bat_fdr] = fdr_bh(pMap_bat);
    nii.img = single(reshape(pMap_bat_fdr, size(brainMask.img)));
    save_untouch_nii(nii, fullfile(pwd, '..', 'ISCtoolbox', 'bathroom',...
        'wholeBrainPvalsFDR.nii'));

    pMap_kit(isnan(kit_map)) = nan;
    [~, ~, ~, pMap_kit_fdr] = fdr_bh(pMap_kit);
    nii.img = single(reshape(pMap_kit_fdr, size(brainMask.img)));
    save_untouch_nii(nii, fullfile(pwd, '..', 'ISCtoolbox', 'bathroom',...
        'wholeBrainPvalsFDR.nii'));

    pMap_avg(isnan(avg_map)) = nan;
    [~, ~, ~, pMap_avg_fdr] = fdr_bh(pMap_avg);
    nii.img = single(reshape(pMap_avg_fdr, size(brainMask.img)));
    save_untouch_nii(nii, fullfile(pwd, '..', 'ISCtoolbox', ...
        'wholeBrainPvalsAverageFDR.nii'));

    % tresholded average map based on fdr correction
    treshold_map = nan(size(avg_map));
    treshold_map(pMap_avg_fdr<0.05) = avg_map(pMap_avg_fdr<0.05);
    nii.img = single(reshape(treshold_map, size(brainMask.img)));
    save_untouch_nii(nii, fullfile(pwd, '..', 'ISCtoolbox', ...
        'wholeBrainAverageTresholded.nii'));

    % treshold on p < 0.001 and 30 voxel cluster
    treshold_p_3d = double(reshape((pMap_avg<0.001), size(brainMask.img)));
    [labels, num] = spm_bwlabel(treshold_p_3d, 18);

    treshold_p_1d = zeros(size(avg_map));
    for c = 1:num
        voxels_in_cluster = find(labels == c);
        if numel(voxels_in_cluster) >= 30
            treshold_p_1d(voxels_in_cluster) = 1;
        end
    end

    treshold_p_1d(treshold_p_1d == 1) = avg_map(treshold_p_1d == 1);
    nii.img = single(reshape(treshold_p_1d, size(brainMask.img)));
    save_untouch_nii(nii, fullfile(pwd, '..', 'ISCtoolbox', ...
        'wholeBrainAverageTresholded_p001_cluster30.nii'));


    % store max values in d
    d.maxR_bat = maxR_bat;
    disp(['Bathroom: Treshold form max r vals: ', num2str(prctile(maxR_bat, 95))])
    d.maxR_kit = maxR_kit;
    disp(['Kitchen: Treshold form max r vals: ', num2str(prctile(maxR_kit, 95))])
    d.maxR_avg = maxR_avg;
    disp(['Category average: Treshold form max r vals: ', num2str(prctile(maxR_avg, 95))])


end

