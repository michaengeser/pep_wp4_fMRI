function results = plot_voxel_steps_per_ROI(cfg, d)

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
if ~isfield(cfg, 'rois_of_interest'); cfg.rois_of_interest = {'LOC', 'PPA', 'TOS', 'LPFC'}; end
if ~isfield(cfg, 'numVoxels'); cfg.numVoxels = [1,2,5,10,20,50,100,200,500,1000,5000,inf]; end % inf = all voxels
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

nROIs = length(cfg.rois_of_interest);
nVoxNums = length(cfg.numVoxels);
allR = nan(nVoxNums, nROIs);
allP = allR;


for iRoi=1:nROIs
    mask_label_short = cfg.rois_of_interest{iRoi};

    % ---------------- Load ISC data ----------------

    fn_bat = fullfile(cfg.outputPath, 'group_level', sprintf('ISC_%s_%s_voxel_steps.mat', mask_label_short, 'bathroom'));
    fn_kit = fullfile(cfg.outputPath, 'group_level', sprintf('ISC_%s_%s_voxel_steps.mat', mask_label_short, 'kitchen'));

    Sbat = load(fn_bat,'meanISC','subPairs');
    Skit = load(fn_kit,'meanISC','subPairs');

    % ---------------- Observed correlations ----------------
    fprintf('Computing observed parcel correlations...\n');
    r_bath = nan(1,nVoxNums);
    r_kit  = nan(1,nVoxNums);

    for iVox = 1:nVoxNums
        v_bat = Sbat.meanISC(:,iVox);
        v_kit = Skit.meanISC(:,iVox);
        if sum(isnan(v_bat)) > 1 || sum(isnan(v_kit)) > 1
            warning(['Parcel number: ' num2str(iVox) ' has NaNs'])
        end
        r_bath(iVox) = corr(v_bat, bat_resid, 'row', 'pairwise', 'type', partial_correlation_type);
        r_kit(iVox)  = corr(v_kit, kit_resid, 'row', 'pairwise', 'type', partial_correlation_type);
    end

    r_avg = (r_bath + r_kit)/2;

    % ---------------- Permutation test ----------------
    fprintf('Running permutations (%d)...\n', nPerm);
    perm_values = zeros(nPerm, nVoxNums, 'single');

    parfor perm_i = 1:nPerm
        V_bat = Sbat.meanISC(random_RDM_vecs{perm_i}, :); % permuted ISC
        V_kit = Skit.meanISC(random_RDM_vecs{perm_i}, :);

        r_bat_perm = nan(1,nVoxNums);
        r_kit_perm = nan(1,nVoxNums);

        for iVox =1:nVoxNums
            vb = V_bat(:,iVox );
            vk = V_kit(:,iVox );
            r_bat_perm(iVox ) = corr(vb, bat_resid, 'row', 'pairwise', 'type', partial_correlation_type);
            r_kit_perm(iVox ) = corr(vk, kit_resid, 'row', 'pairwise', 'type', partial_correlation_type);
        end

        perm_values(perm_i,:) = single((r_bat_perm + r_kit_perm)/2);
    end

    % ---------------- p-values ----------------
    fprintf('Computing empirical p-values...\n');
    p_uncorr_perm = nan(1,nVoxNums);
    for iVox =1:nVoxNums
        obs = r_avg(iVox);
        null = perm_values(:,iVox);
        p_uncorr_perm(iVox) = (sum(null >= obs)+1)/(nPerm+1); % right tailed
    end

    % get p values
    p_uncorr = p_uncorr_perm;
    [~,~,~,p_fdr] = fdr_bh(p_uncorr);

    % ---------------- Return & save ----------------
    results(iRoi).all.r_bathroom = r_bath;
    results(iRoi).all.r_kitchen  = r_kit;
    results(iRoi).all.r_avg      = r_avg;
    results(iRoi).all.p_uncorr   = p_uncorr;
    results(iRoi).all.p_fdr      = p_fdr;
    results(iRoi).ROI            = mask_label_short;
    if cfg.save_perms, results(iRoi).all.perm_values = perm_values; end

    allR(:, iRoi) = r_avg;
    allP(:, iRoi) = p_fdr;
end

% Save results
outfn = fullfile(cfg.outputPath, 'group_level', 'searchlight_parcel_ISC_results.all.mat');
save(outfn,'results','-v7.3');
fprintf('Saved results to %s\n', outfn);

% plot results
    clrs = [
        1.00, 0.00, 1.00;  % LOC
        0.73, 0.02, 0.73;  % PPA
        0.47, 0.03, 0.45;  % TOS
        0.24, 0.04, 0.30;  % LPFC
        ];

figure;
hold on

allP = allP';
allR = allR';

for i = 1:height(allP)

    % plot line
    allLinePlots(i) = plot(allR(i, :)', 'Color', clrs(i, :),  'LineWidth',2);
    for ii = 1:width(allP)

        % add significance markers
        if allP(i, ii) < 0.05
            plot(ii, allR(i, ii), 'o', 'MarkerFaceColor', clrs(i, :), 'MarkerEdgeColor', clrs(i, :))
        end

    end
end

% get x ticks
xticks(1:length(cfg.numVoxels))
xticklabels(string(cfg.numVoxels))
if 1 < length(cfg.numVoxels)
    xlim([1, length(cfg.numVoxels)])
end
xtickangle(45)
xlabel('Number of voxels')

    legendLables = cfg.rois_of_interest;
    legend(allLinePlots, legendLables, 'Location','eastoutside')
    ylim([-0.1, 0.2])


set(gca, 'LineWidth', 2, 'FontName', cfg.FontName, 'FontSize', cfg.FontSize, 'FontWeight', 'bold')
ax = gca;
ax.Box = 'off';

end
