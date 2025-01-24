function cfg = defaultCfg(cfg)

%% set default values for cfg

% exp parameters
if ~isfield(cfg, 'subNums'); cfg.subNums = 2:5;end
if ~isfield(cfg, 'n'); cfg.n = length(cfg.subNums);end
if ~isfield(cfg, 'categories'); cfg.categories = {'bathroom', 'kitchen'};end
if ~isfield(cfg, 'exp_num'); cfg.exp_num = 1;end
if ~isfield(cfg, 'nTrials'); cfg.nTrials = 100;end
if ~isfield(cfg, 'nTargets'); cfg.nTargets = 16;end

% frmi parameter
if ~isfield(cfg, 'nRuns'); cfg.nRuns = 10; end
if ~isfield(cfg, 'tr'); cfg.tr = 1.85;end
if ~isfield(cfg, 'nVols'); cfg.nVols = 152;end
if ~isfield(cfg, 'map'); cfg.map = 't';end
if ~isfield(cfg, 'rois'); cfg.rois = {'wV1.nii', 'wV2.nii',...
        'wLOC.nii', 'wPPA.nii',...
        'wTOS.nii', 'wRSC.nii'};
end

% analsyis paramters
if ~isfield(cfg, 'dnn'); cfg.dnn = 'vgg16_imagenet';end
if ~isfield(cfg, 'functionPath'); cfg.functionPath = fullfile(pwd,'utilities');end
if ~isfield(cfg, 'analysis_names'); cfg.analysis_names = {'typical', 'control'};end
if ~isfield(cfg, 'colormaps'); cfg.colormaps = load(fullfile(cfg.functionPath, 'white_zero_colormap.mat'));end
if ~isfield(cfg, 'plotting'); cfg.plotting = true; end
if ~isfield(cfg, 'saving'); cfg.saving = true; end




end
