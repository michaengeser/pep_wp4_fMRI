function allSplitCorr = iterateSplits(allSplitCorr, splits, iRoi, category, d, cfg)

% preparation
if ~isfield(cfg, 'is_permutation_test'); cfg.is_permutation_test = false;end

for iSplit = 1:height(splits)

    % get first half
    currentRDM = {d.([category,'_RDM']).timecourseRDM(iRoi).runwiseRDM(splits(iSplit,:)).RDM};
    currentRDM = cell2mat(currentRDM);
    currentRDM = reshape(currentRDM, cfg.n, cfg.n, []);
    mean1stHalf = mean(currentRDM, 3, 'omitnan');
    if cfg.is_permutation_test
        mean1stHalf = mean1stHalf(cfg.random_seq, cfg.random_seq);
    end
    mean1stHalf(eye(cfg.n) == 1) = 0;
    mean1stHalf = squareform(mean1stHalf);

    % get second half
    otherRuns = setdiff(1:cfg.nRuns/2, splits(iSplit,:));
    currentRDM = {d.([category,'_RDM']).timecourseRDM(iRoi).runwiseRDM(otherRuns).RDM};
    currentRDM = cell2mat(currentRDM);
    currentRDM = reshape(currentRDM, cfg.n, cfg.n, []);
    mean2ndHalf = mean(currentRDM, 3, 'omitnan');
    mean2ndHalf(eye(cfg.n) == 1) = 0;
    mean2ndHalf = squareform(mean2ndHalf);

    % get correaltiom
    r = corr([mean1stHalf', mean2ndHalf'], 'Type', 'Pearson');
    allSplitCorr(iSplit) = r(1,2);

end
end