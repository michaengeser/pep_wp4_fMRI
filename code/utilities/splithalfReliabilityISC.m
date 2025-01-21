function d = splithalfReliabilityISC(d, cfg)



% evaluate input
if ~isfield(cfg, 'plotting'); cfg.plotting = true; end
if ~isfield(cfg, 'saving'); cfg.saving = false; end
if ~isfield(cfg, 'nRuns'); cfg.nRuns = 10; end

% get combination to split the 10 runs into 2 halfs
splits = nchoosek(1: cfg.nRuns, cfg.nRuns/2);
splits = splits(1:height(splits)/2, :);

% loop through categories
for category = cfg.categories
    category = char(category);
    Category = strcat(upper(category(1)),lower(category(2:end))); % capitalize first letter

    % loop through ROIs
    allRoiCorr = nan(1, numel(cfg.rois));
    for iRoi = 1:numel(cfg.rois)
        roi = char(cfg.rois{iRoi});
        mask_label_short = split(roi, '.');
        mask_label_short = mask_label_short{1};

        % looop through splits
            allSplitCorr = nan(1, height(splits));
            for iSplit = 1:height(splits)

                % get first half
                currentRDM = {d.([category,'_RDM']).timecourseRDM(iRoi).runwiseRDM(splits(iSplit,:)).RDM};
                currentRDM = cell2mat(currentRDM);
                currentRDM = reshape(currentRDM, cfg.n, cfg.n, []);
                mean1stHalf = mean(currentRDM, 3);
                mean1stHalf(eye(cfg.n) == 1) = 0;
                mean1stHalf = squareform(mean1stHalf);

                % get second half
                otherRuns = setdiff(1:cfg.nRuns, splits(iSplit,:));
                currentRDM = {d.([category,'_RDM']).timecourseRDM(iRoi).runwiseRDM(otherRuns).RDM};
                currentRDM = cell2mat(currentRDM);
                currentRDM = reshape(currentRDM, cfg.n, cfg.n, []);
                mean2ndHalf = mean(currentRDM, 3);
                mean2ndHalf(eye(cfg.n) == 1) = 0;
                mean2ndHalf = squareform(mean2ndHalf);

                % get correaltiom
                r = corr([mean1stHalf', mean2ndHalf'], 'Type', 'Spearman');
                allSplitCorr(iSplit) = r(1,2);

            end 

            % take mean across splits
            allRoiCorr(iRoi) = mean(allSplitCorr);
            d.splitHalfISC.(Category).(mask_label_short).allSplitCorr = allSplitCorr;

    end % roi loop

    % write to struct
    d.splitHalfISC.(Category).(mask_label_short).meanCorr = allRoiCorr;


    %% plotting

    if cfg.plotting

        % bar plot
        figure;
        bar(allRoiCorr)
        ylabel('Split Half Correlation')
        xlabel('ROI')
        xticklabels(cfg.rois)
        ylim([min(allRoiCorr) - 0.01, max(allRoiCorr) + 0.01])
        title(['Split-half Reliability of ISC across Runs - ', category])
    end 
end % category loop


end