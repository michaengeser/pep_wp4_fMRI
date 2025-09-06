function d = splithalfReliabilityISC(d, cfg)

% evaluate input
if ~isfield(cfg, 'plotting'); cfg.plotting = true; end
if ~isfield(cfg, 'saving'); cfg.saving = false; end
if ~isfield(cfg, 'nRuns'); cfg.nRuns = 12; end

% get combination to split the 10 runs into 2 halfs
splits = nchoosek(1: cfg.nRuns/2, cfg.nRuns/4);
splits = splits(1:height(splits)/2, :);

% loop through categories
for category = cfg.categories
    category = char(category);
    Category = strcat(upper(category(1)),lower(category(2:end))); % capitalize first letter

    % run time control
    disp(['Evaluate ', category])

    % loop through ROIs
    allRoiCorr = nan(1, numel(cfg.rois));
    for iRoi = 1:numel(cfg.rois)
        roi = char(cfg.rois{iRoi});
        mask_label_short = split(roi, '.');
        mask_label_short = mask_label_short{1};

        % run time control
        disp(['   - ', mask_label_short])

        % looop through splits
        cfg.is_permutation_test = false;
        allSplitCorr = nan(1, height(splits));
        allSplitCorr = iterateSplits(allSplitCorr, splits, iRoi, category, d, cfg);

        % take mean across splits
        allRoiCorr(iRoi) = mean(allSplitCorr, 'omitnan');
        d.splitHalfISC.(Category).(mask_label_short).allSplitCorr = allSplitCorr;

        % run permutation test
        if cfg.permutation_test
            rng('default')
            rng(1)

            for perm = 1:cfg.n_permutations

                % get randiom sequence and performe correaltion
                cfg.random_seq = randperm(cfg.n);
                cfg.is_permutation_test = true;
                allSplitCorr = nan(1, height(splits));
                allSplitCorr = iterateSplits(allSplitCorr, splits, iRoi, category, d, cfg);
                perm_r_vals(perm) = mean(allSplitCorr, 'omitnan');
            end
        end

        permRes(iRoi).(category)= perm_r_vals;

    end % roi loop

    % write to struct
    d.splitHalfISC.(Category).meanCorr = allRoiCorr;

    %% plotting

%     if cfg.plotting
% 
%         % bar plot
%         figure;
%         bar(allRoiCorr)
%         ylabel('Split Half Correlation')
%         xlabel('ROI')
%         xticklabels(cfg.rois)
% 
%         if (min(allRoiCorr) - 0.01) > 0
%             yMin = 0;
%         else
%             yMin = min(allRoiCorr);
%         end 
% 
%         ylim([yMin, max(allRoiCorr) + 0.01])
%         title(['Split-half Reliability of ISC across Runs - ', category])
%     end
end % category loop


if cfg.plotting

% get mean over categories
        mean_cor = (d.splitHalfISC.Bathroom.meanCorr + d.splitHalfISC.Kitchen.meanCorr)/2;

        if cfg.permutation_test
            % get p values and confidence intervals from permutaiotn
            mean_perm_cors = nan(numel(cfg.rois), cfg.n_permutations);
            mean_perm_p = nan(1, numel(cfg.rois));
            mean_perm_ci_upper = mean_perm_p;
            mean_perm_ci_lower = mean_perm_p;
            for iRoi = 1:numel(cfg.rois)
                mean_perm_cors(iRoi, :) = (permRes(iRoi).bathroom +...
                    permRes(iRoi).kitchen)/2;
                mean_perm_p(iRoi) = sum(mean_perm_cors(iRoi, :) >= mean_cor(iRoi))/ cfg.n_permutations;
                mean_perm_ci_upper(iRoi) = mean_cor(iRoi) - prctile(mean_perm_cors(iRoi, :), 5);
                mean_perm_ci_lower(iRoi) = mean_cor(iRoi) - prctile(mean_perm_cors(iRoi, :), 95);
            end

            % get asterisks
            asterisks = pval2asterisks(mean_perm_p, 'fdr');
        end

    % bar plot
    figure;
    hold on
    clrs = bone(numel(cfg.rois)+4);
    clrs = clrs(5:end, :);
    clrs(:, 1:2) = clrs(:, 1:2).*clrs(:, 1:2);
    for iRoi = 1:numel(cfg.rois)
        bar(iRoi, mean_cor(iRoi), 'FaceColor', clrs(iRoi, :))
    end

    % add significance
    if cfg.permutation_test
        for iRoi = 1:numel(cfg.rois)
            errorbar(iRoi, mean_cor(iRoi), ...
                mean_cor(iRoi) - mean_perm_ci_lower(iRoi), ...
                mean_perm_ci_upper(iRoi) - mean_cor(iRoi),...
                'k', 'LineWidth', 1);  % Error bars

            text(iRoi, mean_perm_ci_upper(iRoi) + 0.03, asterisks{iRoi},...
                'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
        end
    end

    % aesthetics
    ylabel(['Pearson correaltion', newline])
    xlabel('ROI')
    xticks(1:numel(cfg.rois))
    xticklabels(cfg.rois)
    xtickangle(45);
    set(gca, 'LineWidth', 1, 'FontName', cfg.FontName, 'FontSize', cfg.FontSize, 'FontWeight', 'bold')
    ax = gca;
    ax.Box = 'off';

    if (min(allRoiCorr) - 0.01) > 0
        yMin = 0;
    else
        yMin = min(allRoiCorr);
    end

    ylim([yMin, max(mean_perm_ci_upper) + 0.1])
    title('Split-half Reliability of ISC across Runs')
end

end