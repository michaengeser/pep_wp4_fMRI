function cfg = defaultCfg(cfg)

%% set default values for cfg

% exp parameters
if ~isfield(cfg, 'subNums'); cfg.subNums = 101:102;end
cfg.n = length(cfg.subNums);
if ~isfield(cfg, 'categories'); cfg.categories = {'bathroom', 'kitchen'};end
if ~isfield(cfg, 'exp_num'); cfg.exp_num = 1;end
if ~isfield(cfg, 'nTrials'); cfg.nTrials = 100;end
if ~isfield(cfg, 'nTargets'); cfg.nTargets = 0;end

% frmi parameter
if ~isfield(cfg, 'nRuns'); cfg.nRuns = 12; end
if ~isfield(cfg, 'stimDur'); cfg.stimDur = 0.25; end
if ~isfield(cfg, 'iti'); cfg.iti = 2; end
if ~isfield(cfg, 'tr'); cfg.tr = 1.85;end
if ~isfield(cfg, 'nVols'); cfg.nVols = 188;end
if ~isfield(cfg, 'map'); cfg.map = 'b';end
if ~isfield(cfg, 'rois'); cfg.rois = {'wV1.nii','wLOC.nii',...
         'wPPA.nii', 'wTOS.nii'};
end

% analsyis paramters
if ~isfield(cfg, 'dnn'); cfg.dnn = 'vgg16_imagenet';end
if ~isfield(cfg, 'functionPath'); cfg.functionPath = fullfile(pwd,'utilities');end
if ~isfield(cfg, 'analysis_names'); cfg.analysis_names = {'typical', 'control'};end
if ~isfield(cfg, 'colormaps'); cfg.colormaps = load(fullfile(cfg.functionPath, 'white_zero_colormap.mat'));end
if ~isfield(cfg, 'plotting'); cfg.plotting = true; end
if ~isfield(cfg, 'saving'); cfg.saving = true; end


% define paths
cfg.sourcedataPath = fullfile(pwd, '..','sourcedata');
cfg.outputPath = fullfile(pwd, '..','derivatives');
cfg.locPath = fullfile(pwd, '..', 'localizer');
cfg.spmPath = (fullfile(pwd, '..', '..', '..', 'MATLAB', 'spm12', 'tpm'));


% other standard configurations
if ~isfield(cfg, 'FontName'); cfg.FontName = 'Helvetica'; end
if ~isfield(cfg, 'FontSize'); cfg.FontSize = 20; end

% load colormap
load('utilities/colormaps.mat');
cfg.colormaps = colormaps;



end
