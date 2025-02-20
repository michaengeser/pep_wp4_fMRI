
voxel_nums = [1,2,3,4,5,10,50,100,150,200];



% preparation
cfg.plotting = false;
cfg.rois = {'wTOS.nii'};
cfg.regressOutMean = true;
cfg.correlation_type = 'Pearson';
cfg.partial_cor = true;
cfg.plot_rdm =false;
cfg.predictor_RDMs = {'typical_late', 'control_late', 'photos_late'};
cfg.RDM_to_partial_out = cfg.predictor_RDMs;
cfg.rois_of_interest = {'TOS'};
cfg.dnns = {'vgg16_imagenet'};
cfg.permutation_test = false;
cfg.bootstrapping = false;
cfg.ISC_type = 'timecourseRDM';
cfg.partial_correlation_type = 'pearson';
cfg.skipIfExists = false;
cfg.smoothing = true;
cfg.saveWholeBrain = false;


rValTypKit = struct;
rValTypBat = struct;
rValCtrKit = struct;
rValCtrBat = struct;

for roi = cfg.rois_of_interest
    rValTypKit.(roi{:}) = nan(1, length(voxel_nums));
    rValTypBat.(roi{:}) = nan(1, length(voxel_nums));
    rValCtrKit.(roi{:}) = nan(1, length(voxel_nums));
    rValCtrBat.(roi{:}) = nan(1, length(voxel_nums));
end

for iVoxeln = 1:length(voxel_nums)

    tic
    % get timecourses (average Voxels for each ROI)
    cfg.n_voxels = voxel_nums(iVoxeln);
    funcROIs(cfg)
    averageVoxelsInROI(cfg)

    % get ISC from timecourses
    d = neural_timecourse_intersub_cor(d, cfg);

    % compare to drawings
    d = compare_roi_RDMs_to_predictor_RDMs(d, cfg);

    for roi = cfg.rois_of_interest
        roiname = roi{:};
        rValTypKit.(roiname)(iVoxeln) = d.compare_roi_to_predictor.(roiname).kitchen.r_val(1);
        rValTypBat.(roiname)(iVoxeln) = d.compare_roi_to_predictor.(roiname).bathroom.r_val(1);
        rValCtrKit.(roiname)(iVoxeln) = d.compare_roi_to_predictor.(roiname).kitchen.r_val(2);
        rValCtrBat.(roiname)(iVoxeln) = d.compare_roi_to_predictor.(roiname).bathroom.r_val(2);
    end
    toc
    disp(['Analysis for ', num2str(voxel_nums(iVoxeln)), ' voxel'])

end


% Get all ROI names from the structure fields (assumes all structures have the same fields)
roiNames = fieldnames(rValTypKit);

% Define the number of ROIs
numROIs = numel(roiNames);

% Create a tiled layout
figure;
t = tiledlayout(ceil(sqrt(numROIs)), ceil(sqrt(numROIs)), 'TileSpacing', 'Compact', 'Padding', 'Compact');

% Loop through each ROI and create a plot in a separate tile
for i = 1:numROIs
    roiname = roiNames{i};  % Get current ROI name

    % Get the vectors for the current ROI
    y1 = rValTypKit.(roiname);
    y2 = rValTypBat.(roiname);
    y3 = rValCtrKit.(roiname);
    y4 = rValCtrBat.(roiname);

    % Get the x-axis (iterations)
    x = 1:length(y1);

    % Select the next tile for plotting
    nexttile;

    % Plot the four lines
    hold on;
    plot(x, y1, '-o', 'LineWidth', 1.5, 'DisplayName', 'TypKit', 'Color', [1, 0, 1]);
    plot(x, y2, '-s', 'LineWidth', 1.5, 'DisplayName', 'TypBat', 'Color', [.9, 0, .9]);
    plot(x, y3, '-o', 'LineWidth', 1.5, 'DisplayName', 'CtrKit', 'Color', [.8, .8, .8]);
    plot(x, y4, '-s', 'LineWidth', 1.5, 'DisplayName', 'CtrBat',  'Color', [.7, .7, .7]);

    % Customize the plot
    title(roiname, 'Interpreter', 'none');  % Title as ROI name
    xlabel('n Voxel');
    xticks(x)
    xticklabels(num2str(voxel_nums'))
    ylabel('Partial correaltion r');
    grid off;

    % Add legend only for the first tile to avoid redundancy
    if i == 1
        legend('Location', 'best');
    end
end

% Add a global title for the tiled layout
sgtitle('Analysis with different Voxel selections for Each ROI');
