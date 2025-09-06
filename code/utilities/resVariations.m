function [] = resVariations(d, cfg)

combis = [1, 1, 1;
    1, 1, 0;
    1, 0, 1;
    1, 0, 0;
    0, 1, 1;
    0, 0, 1;
    0, 1, 0;
    0, 0, 0];

for combi = 1:height(combis)

    % get timecourses (average Voxels for each ROI)
    cfg.plotting = false;
    if combis(combi, 1) == 1
        cfg.cutTargets = true;
        target_str = 'cut targets, ';
    else
        cfg.cutTargets = false;
        target_str = 'inlcude targets, ';
    end
    cfg.skipIfExists = true;
    cfg.smoothing = false;
    cfg.saveWholeBrain = false;
    averageVoxelsInROI(cfg)

    % get ISC from timecourses
    cfg.plot_rdm = false;
    cfg.partial_cor = false;
    cfg.correlation_type = 'Pearson';
    if combis(combi, 2) == 1
        cfg.regressOutMean = true;
        regressOutMean_str = 'with mean regression, ';
    else
        cfg.regressOutMean = false;
        regressOutMean_str = 'no mean regression, ';
    end
    if combis(combi, 3) == 1
        cfg.detrend = true;
        detrend_str = 'detrend ';
    else
        cfg.detrend = false;
        detrend_str = 'no detrend ';
    end
    cfg.correlation_type = 'Pearson';
    d = neural_timecourse_intersub_cor(d, cfg);

    % compare to drawings

    %preparation
    cfg.partial_cor = true;
    cfg.plotting = true;
    cfg.plot_rdm = false;
    cfg.predictor_RDMs = {'typical_late', 'control_late', 'photos_late'};
    cfg.RDM_to_partial_out = cfg.predictor_RDMs;
    % cfg.rois_of_interest = {'V1', 'V2', ...
    %     'LOC', 'PPA', 'TOS', 'RSC', 'LPFC',...
    %     'motorCortex', 'S1', 'auditoryCortex'};
    cfg.rois_of_interest = {'V1', 'LOC',...
        'PPA', 'TOS', 'LPFC'};
    cfg.dnns = {'vgg16_imagenet'};
    cfg.plot_type = 'bar';
    cfg.show_single_cate = false;
    cfg.permutation_test = true;
    cfg.n_permutations = 10000;
    cfg.permutation_type = 'row_col_shuffle_ref';
    cfg.ISC_type = 'timecourseRDM';
    cfg.partial_correlation_type = 'Pearson';
    cfg.plotting_predictors = 1;
     cfg.add_legend = false;

    d = compare_roi_RDMs_to_predictor_RDMs(d, cfg);

    title([target_str, regressOutMean_str, detrend_str, newline])
end
end