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
if ~isfield(cfg, 'rois_of_interest'); cfg.rois_of_interest = {'V1', 'LOC', 'PPA', 'TOS', 'LPFC'}; end
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

% preallocate
nParcels = 180;
allR = nan(nParcels, length(cfg.numVoxels));
allP = allR;

for iVox = 1:length(cfg.numVoxels)
    nVox = cfg.numVoxels(iVox);


    % ---------------- Load ISC data ----------------

    fn_bat = fullfile(cfg.outputPath, 'group_level', sprintf('ISC_HCP_%s_%d_voxels.mat', 'bathroom', nVox));
    fn_kit = fullfile(cfg.outputPath, 'group_level', sprintf('ISC_HCP_%s_%d_voxels.mat', 'kitchen', nVox));

    Sbat = load(fn_bat,'meanISC','subPairs');
    Skit = load(fn_kit,'meanISC','subPairs');

    % ---------------- Observed correlations ----------------
    fprintf('Computing observed parcel correlations...\n');
    r_bath = nan(1,nParcels);
    r_kit  = nan(1,nParcels);

    for p=1:nParcels
        v_bat = Sbat.meanISC(:,p);
        v_kit = Skit.meanISC(:,p);
        if sum(isnan(v_bat)) > 1 || sum(isnan(v_kit)) > 1
            warning(['Parcel number: ' num2str(p) ' has NaNs'])
            %             disp(v_bat)
            %             disp(newline)
            %             disp(v_kit)
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
    results(iVox).all.r_bathroom = r_bath;
    results(iVox).all.r_kitchen  = r_kit;
    results(iVox).all.r_avg      = r_avg;
    results(iVox).all.p_uncorr   = p_uncorr;
    results(iVox).all.p_fdr      = p_fdr;
    results(iVox).nVox           = nVox;
    if cfg.save_perms, results(iVox).all.perm_values = perm_values; end

    allR(:, iVox) = r_avg';
    allP(:, iVox) = p_fdr';


    %% restrict analysis to parcel subsets

    % Load the atlas labels
    labels = readtable(fullfile(pwd, '..', 'MNI_ROIs', 'HCP_atlas_lables.csv'));
    labels = labels(1:nParcels, :);
    labels.r = r_avg';
    labels.p_fdr_all = p_fdr';
    labels.uncorrected_p = p_uncorr';
    labels.uncorrected_p_norm = p_norm';
    results(iVox).labels_all = labels;

    %
    %         % Requested parcels
    %         parcel_names = {...
    %             'FEF','PEF','55b','8Av','8Ad','9m','8BL','9p','10d','8C',...
    %             '44','45','47l','a47r','6r','IFJa','IFJp','IFSp','IFSa',...
    %             'p9-46v','46','a9-46v','9-46d','9a','10v','a10p','10pp',...
    %             '6a','i6-8','s6-8','p10p','p47r'};
    %
    %         % Match ROI numbers
    %         mask = ismember(labels.roi, parcel_names);
    %         lpfc_rois = labels.num_roi(mask);
    %
    %         % index LPFC
    %         lpfc_rs = r_avg(lpfc_rois);
    %         p_uncorr_lpfc = p_uncorr(lpfc_rois);
    %         [~,~,~,p_fdr_lpfc] = fdr_bh(p_uncorr_lpfc);
    %
    %         % lpfc parcels table
    %         lpfcParcelTable = labels(lpfc_rois, :);
    %         lpfcParcelTable.p_fdr = p_fdr_lpfc';
    %         lpfcParcelTable.r = r_avg(lpfc_rois)';
    %
    %         % significant parcels
    %         sigParcels = lpfc_rois((p_fdr_lpfc < 0.05));
    %         sigParcelsLabels = labels(sigParcels, :);
    %         sigParcelsLabels.p_fdr = p_fdr_lpfc((p_fdr_lpfc < 0.05))';
    %         sigParcelsLabels.r = lpfc_rs((p_fdr_lpfc < 0.05))';
    %         disp(newline)
    %         disp('Significant parcels in lateeral prefrontal cortex')
    %         disp(newline)
    %         disp(sigParcelsLabels)
    %
    %         results.lpfc.p_uncorr = p_uncorr_lpfc;
    %         results.lpfc.p_fdr = p_fdr_lpfc;
    %         results.lpfc.significant_parcels = sigParcels;
    %         results.lpfc.significant_parcels_labels = sigParcelsLabels.roi;
    %
    %
    %         % visual cortex
    %         visual_rois = [...
    %             1,2,3,4,5,6,7,13,16,17,18,19,20,21,22,23,...
    %             137,138,152,153,154,156,157,158,159,160,163];
    %
    %         % index visual cortex
    %         vis_rs = r_avg(visual_rois);
    %         p_uncorr_vis = p_uncorr(visual_rois);
    %         [~,~,~,p_fdr_vis] = fdr_bh(p_uncorr_vis);
    %
    %         % visual parcels table
    %         visParcelTable = labels(visual_rois, :);
    %         visParcelTable.p_fdr = p_fdr_vis';
    %         visParcelTable.r = r_avg(visual_rois)';
    %
    %         % significant parcels
    %         sigParcels = visual_rois((p_fdr_vis < 0.05));
    %         sigParcelsLabels = labels(sigParcels, :);
    %         sigParcelsLabels.p_fdr = p_fdr_vis((p_fdr_vis < 0.05))';
    %         sigParcelsLabels.r = vis_rs((p_fdr_vis < 0.05))';
    %         disp(newline)
    %         disp('Significant parcels in visual cortex')
    %         disp(newline)
    %         disp(sigParcelsLabels)
    %
    %         results.visual.p_uncorr = p_uncorr_vis;
    %         results.visual.p_fdr = p_fdr_vis;
    %         results.visual.significant_parcels = sigParcels;
    %         results.visual.significant_parcels_labels = sigParcelsLabels.roi;


end

% Save results
outfn = fullfile(cfg.outputPath, 'group_level', 'searchlight_parcel_ISC_results.all.mat');
save(outfn,'results','-v7.3');
fprintf('Saved results to %s\n', outfn);

% plot results

%clrs = jet(height(allP));
clrMap = cfg.colormaps.gist_ncar(1:end-10, :);
clrs = clrMap(round((1:nParcels) * (height(clrMap)/nParcels)), :);


figure;
hold on

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
xlabel('Number of voxels')
title('All parcels')
set(gca, 'LineWidth', 2, 'FontName', cfg.FontName, 'FontSize', cfg.FontSize, 'FontWeight', 'bold')
ax = gca;
ax.Box = 'off';

% get all parcels with at least one sigificant test

% add labels

sigPoints = allP < 0.05;
sigParcels = find(mean(sigPoints, 2) > 0);
sigR = allR(sigParcels, :);
sigP = allP(sigParcels, :);

% plot results
figure;
hold on

for i = 1:height(sigP)

    % plot line
    linePlots(i) = plot(sigR(i, :)', 'Color', clrs(sigParcels(i), :));
    for ii = 1:width(sigP)

        % add significance markers
        if sigP(i, ii) < 0.05
            plot(ii, sigR(i, ii), 'o',...
                'MarkerFaceColor', clrs(sigParcels(i), :), 'MarkerEdgeColor', clrs(sigParcels(i), :))
        end

    end
end

% configure axis
xticks(1:length(cfg.numVoxels))
xticklabels(string(cfg.numVoxels))
xlim([1, length(cfg.numVoxels)])
xlabel('Number of voxels')
xtickangle(45)
ylabel(['Partial correlation [r]', newline]);
ylim([-0.1, 0.3])
yline(0, 'LineWidth', 2, 'Color', 'k');

% add labels
labels = readtable(fullfile(pwd, '..', 'MNI_ROIs', 'HCP_atlas_lables.csv'));
legendLables = labels.community(sigParcels);
legend(linePlots, legendLables, 'Location','eastoutside')


title('Parcels with at least one significant test')


%% plot all parcel that are significant at 200 voxels
sig200Points = allP(:, cfg.numVoxels == 200) < 0.05;
sig200Parcels = find(sig200Points > 0);
sig200R = allR(sig200Parcels, :);
sig200P = allP(sig200Parcels, :);
sig200Plabels = labels(sig200Parcels, :);

% make new Order for plotting
newOrder = [3,7,8,9,5,6,4,10,11,13,1,12,2]; % hard coded
sig200R = sig200R(newOrder, :);
sig200P = sig200P(newOrder, :);
sig200Plabels = sig200Plabels(newOrder, :);

structPos = find(cfg.numVoxels == 200);
results(structPos).sig200Parcels = sig200Parcels(newOrder);
results(structPos).newOrder = newOrder;

% plot results
figure;
hold on

% get colors
load(fullfile(pwd, 'utilities', 'parcelColors.mat'))

for i = 1:height(sig200P)
    % plot line
    linePlots(i) = plot(sig200R(i, :)', 'Color', parcelColors{i}, 'LineWidth',2);
    legendLabel{i} = [sig200Plabels.community{i}, ' (', sig200Plabels.roi{i}, ')'];
    for ii = 1:width(sig200P)
        % add significance markers
        if sig200P(i, ii) < 0.05
            plot(ii, sig200R(i, ii), 'o',...
                'MarkerFaceColor', [parcelColors{i}], 'MarkerEdgeColor', parcelColors{i})
        end
    end
end

% get x ticks
xticks(1:length(cfg.numVoxels))
xticklabels(string(cfg.numVoxels))
xtickangle(45)
xlim([1, length(cfg.numVoxels)])
xlabel('Number of voxels')
yline(0, 'LineWidth', 2, 'Color', 'k');

% add labels
labels = readtable(fullfile(pwd, '..', 'MNI_ROIs', 'HCP_atlas_lables.csv'));
legend(linePlots, legendLabel, 'Location','eastoutside')


set(gca, 'LineWidth', 2, 'FontName', cfg.FontName, 'FontSize', cfg.FontSize, 'FontWeight', 'bold')
ax = gca;
ax.Box = 'off';

end
