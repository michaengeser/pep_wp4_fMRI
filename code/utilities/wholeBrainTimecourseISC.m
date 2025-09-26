function d = wholeBrainTimecourseISC(d, cfg)
% evaluate input
if ~isfield(cfg, 'plotting'); cfg.plotting = true; end
if ~isfield(cfg, 'saving'); cfg.saving = false; end
if ~isfield(cfg, 'nRuns'); cfg.nRuns = 12; end
if ~isfield(cfg, 'regressOutMean'); cfg.regressOutMean = true; end
if ~isfield(cfg, 'brainMask'); cfg.brainMask = 'group_mask_common'; end
if ~isfield(cfg, 'predictor_RDMs'); cfg.predictor_RDMs = {'typical_late', 'control_late', 'photos_late'}; end
if ~isfield(cfg, 'RDM_to_partial_out'); cfg.RDM_to_partial_out = cfg.predictor_RDMs; end
if ~isfield(cfg, 'correlation_type'); cfg.correlation_type = 'pearson';end
if ~isfield(cfg, 'permutation_test'); cfg.permutation_test = true;end
if ~isfield(cfg, 'n_permutations'); cfg.n_permutations = 10000;end
if ~isfield(cfg, 'partial_correlation_type'); cfg.partial_correlation_type = 'Pearson';end
if ~isfield(cfg, 'dnns'); cfg.dnns = {cfg.dnn};end
if ~isfield(cfg, 'p_thresh'); cfg.p_thresh = 0.001;end
if ~isfield(cfg, 'smoothKernel'); cfg.smoothKernel = 6;end

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
cfg.RDM_to_partial_out = cfg.predictor_RDMs;

% convert to local configurations
if cfg.plotting; plotting = true; else plotting = false; end
if cfg.saving; saving = true; else saving = false; end
cfg.plotting = false;
cfg.saving = false;

% get brain mask
if contains(cfg.brainMask, 'group_mask')
    brainMaskPath = fullfile(pwd, '..', 'derivatives', 'group_level', [cfg.brainMask, '.nii']);
    if ~exist(brainMaskPath, 'file')
        makeBrainMask(cfg)
    end
    brainMask = load_untouch_nii(brainMaskPath);
    brainMaskImg = brainMask.img;
elseif strcmp(cfg.brainMask, 'cortexHPC')
    brainMaskPath = fullfile(pwd, '..', 'MNI_ROIs', 'wHPC_atlas.nii');
    brainMask = load_untouch_nii(brainMaskPath);
    brainMaskImg = 0 < brainMask.img;
else 
    brainMaskPath = fullfile(pwd, '..', 'MNI_ROIs', [cfg.brainMask, '.nii']);
    brainMask = load_untouch_nii(brainMaskPath);
    brainMaskImg = 0 < brainMask.img;
end


% load brain mask
% loop through categories
for category = cfg.categories
    category = char(category);
    Category = strcat(upper(category(1)),lower(category(2:end))); % capitalize first letter

    % ISC path
    % ISC path
    if cfg.smoothKernel == 12
        ISCpath = fullfile(pwd, '..', 'ISCtoolbox', [category, '2_s12'], 'results', 'memMaps.mat');
    elseif cfg.smoothKernel == 6
        ISCpath = fullfile(pwd, '..', 'ISCtoolbox', [category, ''], 'results', 'memMaps.mat');
    end
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
isc_values_bat = reshapedxyzc.bathroom';
isc_values_kit = reshapedxyzc.kitchen';

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

if isempty(gcp('nocreate'))
    parpool(5);
end
notNaNvec = true(size(isc_values_bat, 2), 1);
parfor iVoxel = 1:size(isc_values_bat, 2)

    % make nan if mean ISC = 0 for this voxel
    if mean(isc_values_bat(:, iVoxel)) == 0 || ...
            mean(isc_values_kit(:, iVoxel)) == 0 || ...
            brainMaskImg(iVoxel) ~= 1
        isc_values_bat(:, iVoxel) = nan;
        isc_values_kit(:, iVoxel) = nan;
        notNaNvec(iVoxel) = false;
    end
end

%% get correlation maps and save them

% write nifi file
nii = brainMask;
nii.hdr.dime.datatype = 16;
nii.hdr.dime.bitpix = 32;
% reset scaling
nii.hdr.dime.scl_slope = 1;
nii.hdr.dime.scl_inter = 0;

% bathroom
bat_map = nan(size(isc_values_bat, 2), 1);
bat_map(notNaNvec) = corr(isc_values_bat(:, notNaNvec), bat_resid, 'Tail', 'right', 'Type', partial_correlation_type);
nii.img = single(reshape(bat_map, size(brainMask.img)));
save_untouch_nii(nii, fullfile(pwd, '..', 'ISCtoolbox', 'bathroom',...
    'wholeBrainCorMap.nii'));

% kitchen
kit_map = nan(size(isc_values_kit, 2), 1);
kit_map(notNaNvec) = corr(isc_values_kit(:, notNaNvec), kit_resid, 'Tail', 'right', 'Type', partial_correlation_type);
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

    % prepare values for parallel processing
    random_RDM_vecs = cfg.random_RDM_vecs;
    n_permutations = cfg.n_permutations;
    isc_values_bat_parfor = isc_values_bat(:, notNaNvec);
    isc_values_kit_parfor = isc_values_kit(:, notNaNvec);
    all_perm_avg_r = single(nan(size(isc_values_kit_parfor, 2), n_permutations));

    parfor perm = 1:n_permutations

        % bathroom
        perm_bat_map = corr(isc_values_bat_parfor(random_RDM_vecs{perm}, :),...
            bat_resid, 'Tail', 'right', 'Type', partial_correlation_type);
        perm_bat_map(isnan(perm_bat_map)) = 0;

        % kitchen
        perm_kit_map = corr(isc_values_kit_parfor(random_RDM_vecs{perm}, :),...
            kit_resid, 'Tail', 'right', 'Type', partial_correlation_type);
        perm_kit_map(isnan(perm_kit_map)) = 0;

        % average
        all_perm_avg_r(:, perm) = (perm_bat_map + perm_kit_map)/2;


        % progress report
        if mod(perm, 1000) == 0
            disp([char(datetime), ' - ', num2str(perm/n_permutations*100), '% of permutations'])
        end

    end

    % get treshold for significance
    thr_vals = nan(size(avg_map));
    thr_vals(notNaNvec) = quantile(all_perm_avg_r, 1-cfg.p_thresh, 2);
    above_thr = avg_map > thr_vals;

    % get p values
    p_vals = nan(size(avg_map));
    p_vals(notNaNvec) = mean(avg_map(notNaNvec) <= all_perm_avg_r, 2, 'omitnan');
    p_vals(isnan(avg_map)) = nan;

    % write nifti of pvals
    nii.img = single(reshape(p_vals, size(brainMask.img)));
    save_untouch_nii(nii, fullfile(pwd, '..', 'ISCtoolbox', ...
        'wholeBrainPvalsAverage.nii'));

    % get p values FDR-corrected and write nifti
    p_vals_fdr = nan(size(p_vals));
    [~, ~, ~, p_vals_fdr(~isnan(p_vals))] = fdr_bh(p_vals(~isnan(p_vals)));
    nii.img = single(reshape(p_vals_fdr, size(brainMask.img)));
    save_untouch_nii(nii, fullfile(pwd, '..', 'ISCtoolbox', ...
        'wholeBrainPvalsAverageFDR.nii'));

    % tresholded average map based on fdr correction
    treshold_map = nan(size(avg_map));
    treshold_map(p_vals_fdr<0.05) = avg_map(p_vals_fdr<0.05);
    nii.img = single(reshape(treshold_map, size(brainMask.img)));
    save_untouch_nii(nii, fullfile(pwd, '..', 'ISCtoolbox', ...
        'wholeBrainAverageTresholdedFDR.nii'));

    % tresholded average map based on uncorrected p values
    treshold_map = nan(size(avg_map));
    treshold_map(p_vals<0.001) = avg_map(p_vals<0.001);
    nii.img = single(reshape(treshold_map, size(brainMask.img)));
    save_untouch_nii(nii, fullfile(pwd, '..', 'ISCtoolbox', ...
        'wholeBrainAverageTresholdedp001.nii'));


    %% make cluster-based tests

    % get treshold for significance
    thr_vals = nan(size(avg_map));
    thr_vals(notNaNvec) = quantile(all_perm_avg_r, 1-cfg.p_thresh, 2);
    thr_vals(thr_vals < 0.1) = 0.1;
    above_thr = avg_map > thr_vals;
    above_thr3D = reshape(above_thr, size(brainMask.img));

    [labels, num] = spm_bwlabel(double(above_thr3D));
    obs_cluster_mass = zeros(1, num);
    num_voxels_in_cluster = zeros(1, num);
    voxels_in_cluster = {};

    for c = 1:num
        voxels_in_cluster{c} = find(labels == c);
        num_voxels_in_cluster(c) = numel(voxels_in_cluster{c});
        obs_cluster_mass(c) = sum(avg_map(voxels_in_cluster{c}));
    end

    % Build null distribution of maximum cluster statistic
    maxClusterSize = zeros(n_permutations, 1);
    maxClusterMass = zeros(n_permutations, 1);

    for p=1:n_permutations
        perm_map = nan(size(avg_map));
        perm_map(notNaNvec) = all_perm_avg_r(:, p);
        perm_mask = perm_map > thr_vals;  % threshold exactly same way

        % if there is clusters over treshold store their stats
        if max(perm_mask) == 1
            perm_mask3D = reshape(perm_mask, size(brainMask.img));
            [labels, num] = spm_bwlabel(double(perm_mask3D));
            obs_cluster_mass_perm = zeros(1, num);
            num_voxels_in_perm_cluster = zeros(1, num);
            for c = 1:num
                voxels_in_perm_cluster = find(labels == c);
                num_voxels_in_perm_cluster(c) = numel(voxels_in_perm_cluster);
                obs_cluster_mass_perm(c) = sum(perm_map(voxels_in_perm_cluster));
            end
            maxClusterSize(p) = max(num_voxels_in_perm_cluster);
            maxClusterMass(p) = max(obs_cluster_mass_perm);
        end
    end

    % cluster p-values (extent-based)
    cluster_pvals_extent = arrayfun(@(s) mean(maxClusterSize >= s), num_voxels_in_cluster);
    sig_clusters_extent = find(cluster_pvals_extent < 0.05);
    cluster_trh = prctile(maxClusterSize, 95);
    disp(['Treshold cluster extent: ', num2str(cluster_trh)])

    % treshold observed correaltions based on clusters (extent)
    treshold_map = zeros(size(avg_map));
    for c = 1:numel(sig_clusters_extent)
        voxels_in_cluster = find(labels == c);
        treshold_map(voxels_in_cluster) = avg_map(voxels_in_cluster);
    end
    nii.img = single(reshape(treshold_map, size(brainMask.img)));
    save_untouch_nii(nii, fullfile(pwd, '..', 'ISCtoolbox', ...
        'wholeBrainAverageTresholdedClusterExtent.nii'));

    % cluster p-values (mass-based)
    cluster_pvals_mass = arrayfun(@(m) mean(maxClusterMass >= m), obs_cluster_mass);
    sig_clusters_mass = find(cluster_pvals_mass < 0.05);
    cluster_trh = prctile(maxClusterMass, 95);
    disp(['Treshold cluster mass: ', num2str(cluster_trh)])

    % treshold observed correaltions based on clusters (extent)
    treshold_map = zeros(size(avg_map));
    for c = 1:numel(sig_clusters_mass)
        voxels_in_cluster = find(labels == c);
        treshold_map(voxels_in_cluster) = avg_map(voxels_in_cluster);
    end
    nii.img = single(reshape(treshold_map, size(brainMask.img)));
    save_untouch_nii(nii, fullfile(pwd, '..', 'ISCtoolbox', ...
        'wholeBrainAverageTresholdedClusterMass.nii'));

    nii.img = single(reshape(thr_vals, size(brainMask.img)));
    save_untouch_nii(nii, fullfile(pwd, '..', 'ISCtoolbox', ...
        'wholeBrainAverageTresholdedTestThrVal.nii'));


end

