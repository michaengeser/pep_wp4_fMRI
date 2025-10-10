function results = searchlight_parcel_ISC(cfg, d)
% SEARCHLIGHT_PARCEL_ISC correlate parcel-wise ISC with predictor RDMs,
% average categories, permutation test, output r, p_uncorr, p_fdr per parcel.

% ---------------- Defaults ----------------
if ~isfield(cfg, 'n'); error('cfg.n (n subjects) required'); end
if ~isfield(cfg, 'subNums'); error('cfg.subNums required'); end
if ~isfield(cfg, 'dnns'); cfg.dnns = {cfg.dnn}; end
if ~isfield(cfg, 'predictor_RDMs'); cfg.predictor_RDMs = {'typical_late','control_late','photos_late'}; end
if ~isfield(cfg, 'n_permutations'); cfg.n_permutations = 10000; end
if ~isfield(cfg, 'outputPath'); cfg.outputPath = fullfile(pwd,'results'); end
if ~isfield(cfg, 'partial_correlation_type'); cfg.partial_correlation_type = 'Pearson';end
if ~exist(cfg.outputPath,'dir'); mkdir(cfg.outputPath); end
if ~isfield(cfg, 'permutation_tail'); cfg.permutation_tail = 'right'; end
if ~isfield(cfg, 'save_perms'); cfg.save_perms = false; end 
if ~isfield(cfg, 'random_RDM_vecs') 
    % generate permutated subjects list rng('default');
    rng(1) % ensure reproducible outcome
    random_seqs = cell(1, cfg.n_permutations);
    RDM_template = squareform(1:nchoosek(cfg.n, 2));
    for i = 1:cfg.n_permutations
        random_seqs{i} = randperm(cfg.n);
        RDM_shuffle = RDM_template(random_seqs{i}, random_seqs{i});
        random_RDM_vecs{i} = squareform(RDM_shuffle);
    end
    cfg.RDM_to_partial_out = cfg.predictor_RDMs;
else
    random_RDM_vecs = cfg.random_RDM_vecs;
end
nPerm = cfg.n_permutations;
partial_correlation_type = cfg.partial_correlation_type;


% ---------------- Predictor residuals ----------------
% Bathroom
RDMs = d.DNN.(cfg.dnn).control.bathroom.subject_mean(1);
labels = {RDMs.name};
[RDMs, cfg.labels] = evaluate_predictor_RDMs(d, RDMs, labels, cfg, 'bathroom');
for iPred = 1:length(cfg.predictor_RDMs)
    RDMs(iPred+1).RDM(eye(cfg.n) == 1) = 0;
    bat_preds(:, iPred) = squareform(RDMs(iPred+1).RDM);
end
X = [bat_preds(:,2:end), ones(size(bat_preds,1),1)];
b = X \ bat_preds(:,1);
bat_resid = bat_preds(:,1) - X*b;

% Kitchen
RDMs = d.DNN.(cfg.dnn).control.kitchen.subject_mean(1);
labels = {RDMs.name};
[RDMs, cfg.labels] = evaluate_predictor_RDMs(d, RDMs, labels, cfg, 'kitchen');
for iPred = 1:length(cfg.predictor_RDMs)
    RDMs(iPred+1).RDM(eye(cfg.n) == 1) = 0;
    kit_preds(:, iPred) = squareform(RDMs(iPred+1).RDM);
end
X = [kit_preds(:,2:end), ones(size(kit_preds,1),1)];
b = X \ kit_preds(:,1);
kit_resid = kit_preds(:,1) - X*b;

% ---------------- Load ISC data ----------------
fn_bat = fullfile(cfg.outputPath, 'group_level', 'ISC_HPC_bathroom.mat');
fn_kit = fullfile(cfg.outputPath, 'group_level', 'ISC_HPC_kitchen.mat');
Sbat = load(fn_bat,'meanISC','subPairs');
Skit = load(fn_kit,'meanISC','subPairs');
nParcels = size(Sbat.meanISC,2);

% ---------------- Observed correlations ----------------
fprintf('Computing observed parcel correlations...\n');
r_bath = nan(1,nParcels);
r_kit  = nan(1,nParcels);

for p=1:nParcels
    v_bat = Sbat.meanISC(:,p);
    v_kit = Skit.meanISC(:,p);
    if sum(isnan(v_bat)) > 1 || sum(isnan(v_kit)) > 1 
        warning(['Parcel number: ' num2str(p) ' has NaNs'])
        disp(v_bat)
        disp(newline)
        disp(v_kit)
    end 
    r_bath(p) = corr(v_bat, bat_resid, 'row', 'pairwise', 'type', partial_correlation_type);
    r_kit(p)  = corr(v_kit, kit_resid, 'row', 'pairwise', 'type', partial_correlation_type);
end

r_avg = (r_bath + r_kit)/2;

% ---------------- Permutation test ----------------
fprintf('Running permutations (%d)...\n', nPerm);
perm_values = zeros(nPerm, nParcels, 'single');

parfor perm_i = 1:nPerm
    V_bat = Sbat.meanISC(random_RDM_vecs{perm_i}, :); % permuted ISC
    V_kit = Skit.meanISC(random_RDM_vecs{perm_i}, :);

    r_bat_perm = nan(1,nParcels);
    r_kit_perm = nan(1,nParcels);

    for p=1:nParcels
        vb = V_bat(:,p);
        vk = V_kit(:,p);
        r_bat_perm(p) = corr(vb, bat_resid, 'row', 'pairwise', 'type', partial_correlation_type);
        r_kit_perm(p) = corr(vk, kit_resid, 'row', 'pairwise', 'type', partial_correlation_type);
    end

    perm_values(perm_i,:) = single((r_bat_perm + r_kit_perm)/2);
end

% ---------------- p-values ----------------
fprintf('Computing empirical p-values...\n');
p_uncorr_perm = nan(1,nParcels);
p_norm = p_uncorr_perm;
for p=1:nParcels
    obs = r_avg(p);
    null = perm_values(:,p);
    p_uncorr_perm(p) = (sum(null >= obs)+1)/(nPerm+1); % right tailed 

     % --- Gaussian parametric p ---
    mu = mean(null);
    sigma = std(null);

    if sigma > 0
        % right-tailed probability under N(mu, sigma^2)
        p_norm(p) = 1-normcdf(obs, mu, sigma); 
    else
        % fallback if no variance in null
        warning('No variance in null')
        p_norm(p) = nan;
    end
end

% get p values
p_uncorr = p_uncorr_perm;
[~,~,~,p_fdr] = fdr_bh(p_uncorr);

% ---------------- Return & save ----------------
results.all.r_bathroom = r_bath;
results.all.r_kitchen  = r_kit;
results.all.r_avg      = r_avg;
results.all.p_uncorr   = p_uncorr;
results.all.p_fdr      = p_fdr;
if cfg.save_perms, results.all.perm_values = perm_values; end

outfn = fullfile(cfg.outputPath, 'group_level', 'searchlight_parcel_ISC_results.all.mat');
save(outfn,'results','-v7.3');
fprintf('Saved results to %s\n', outfn);



%% restrict analysis to parcel subsets


% Load the atlas labels
labels = readtable(fullfile(pwd, '..', 'MNI_ROIs', 'HPC_atlas_lables.csv'));

% Requested parcels
parcel_names = {...
    'FEF','PEF','55b','8Av','8Ad','9m','8BL','9p','10d','8C',...
    '44','45','47l','a47r','6r','IFJa','IFJp','IFSp','IFSa',...
    'p9-46v','46','a9-46v','9-46d','9a','10v','a10p','10pp',...
    '6a','i6-8','s6-8','p10p','p47r'};

% Match ROI numbers
mask = ismember(labels.roi, parcel_names);
lpfc_rois = labels.num_roi(mask);

% index LPFC
lpfc_rs = r_avg(lpfc_rois);
p_uncorr_lpfc = p_uncorr(lpfc_rois);
[~,~,~,p_fdr_lpfc] = fdr_bh(p_uncorr_lpfc);

% significant parcels 
sigParcels = lpfc_rois((p_fdr_lpfc < 0.05));
sigParcelsLabels = labels(sigParcels, :);
sigParcelsLabels.p = p_fdr_lpfc((p_fdr_lpfc < 0.05))';
sigParcelsLabels.r = lpfc_rs((p_fdr_lpfc < 0.05))';
disp(newline)
disp('Significant parcels in lateeral prefrontal cortex')
disp(newline)
disp(sigParcelsLabels)

results.lpfc.p_uncorr = p_uncorr_lpfc;
results.lpfc.p_fdr = p_fdr_lpfc;
results.lpfc.significant_parcels = sigParcels;
results.lpfc.significant_parcels_labels = sigParcelsLabels.roi;


% visual cortex
visual_rois = [...
    1,2,3,4,5,6,7,13,16,17,18,19,20,21,22,23,...
    137,138,152,153,154,156,157,158,159,160,163];

% index visual cortex
vis_rs = r_avg(visual_rois);
p_uncorr_vis = p_uncorr(visual_rois);
[~,~,~,p_fdr_vis] = fdr_bh(p_uncorr_vis);

% significant parcels 
sigParcels = visual_rois((p_fdr_vis < 0.05));
sigParcelsLabels = labels(sigParcels, :);
sigParcelsLabels.p = p_fdr_vis((p_fdr_vis < 0.05))';
sigParcelsLabels.r = vis_rs((p_fdr_vis < 0.05))';
disp(newline)
disp('Significant parcels in visual cortex')
disp(newline)
disp(sigParcelsLabels)

results.visual.p_uncorr = p_uncorr_vis;
results.visual.p_fdr = p_fdr_vis;
results.visual.significant_parcels = sigParcels;
results.visual.significant_parcels_labels = sigParcelsLabels.roi;


end
