function config = get_config(subs, purpose, varargin)
arguments
    subs (1,:) double
    purpose (1,:) char { ...
        mustBeMember( ...
        purpose,{'create_mcf', 'preproc', 'glm', 'rsa', 'fusion'} ...
        ) ...
        }
end

arguments (Repeating)
    varargin
end

% define optional parameters
p = inputParser;
addOptional( ...
    p, 'glm_level', [], @(x) validate_glm_lvl(x) ...
    );
addOptional( ...
    p, 'glm_normalized', [], ...
    @(x) validateattributes(x, {'logical'}, {'nonempty'}) ...
    );
addOptional( ...
    p, 'fusion_type', [], @(x) validate_fusion_type(x) ...
    );
addOptional( ...
    p, 'modality2avg', [], @(x) validate_modality2avg(x) ...
    );
addOptional( ...
    p, 'commonality_corr_type', [], @(x) validate_comm_corr(x) ...
    );
addOptional( ...
    p, 'commonality_searchlight_size', [], ...
    @(x) isempty(x) || @(x) validateattributes(x, {'double'}) ...
    );
addOptional( ...
    p, 'plots_stats', [], ...
    @(x) validateattributes(x, {'logical'}, {'nonempty'}) ...
    );
addOptional( ...
    p, 'save_comm_maps', [], ...
    @(x) validateattributes(x, {'logical'}, {'nonempty'}) ...
    );
parse(p, varargin{:});

% extract optional parameters
glm_level = p.Results.glm_level;
glm_normalized = p.Results.glm_normalized;

% init.
config = struct();

% define paths
config.sourcedataPath = fullfile(pwd, '..','sourcedata');
outputPath = fullfile(pwd, '..','derivatives');
locPath = fullfile(pwd, '..', 'localizer');
behPath = fullfile(config.sourcedataPath, 'beh');
fmriPath = fullfile(config.sourcedataPath, 'fmri/bids');

% define spm12 path
spmPath = (fullfile(pwd, '..', '..', '..', 'MATLAB', 'spm12', 'tpm'));

% set params. depending on the purpose
switch purpose
    case 'create_mcf'
        config.subs = cellstr(string(subs));
        config.behPath = behPath;
        config.locPath = locPath;
    case 'preproc'
        config.subs = cellstr(string(subs));
        config.fmriPath = fmriPath;
        config.outputPath = outputPath;
        config.spmPath = spmPath;
    case 'glm'
        config.subs = cellstr(string(subs));
        config.fmriPath = fmriPath;
        config.behPath = behPath;
        config.outputPath = outputPath;

        % handle optional GLM parameters
        if ~isempty(glm_level)
            config.glm_level = glm_level;
        end
        if ~isempty(glm_normalized)
            config.glm_normalized = glm_normalized;
        end
    case 'rsa'
        config.subs = subs;
        config.fmriPath = fmriPath;
        config.outputPath = outputPath;
        config.modelRDM = get_predictor_RDM(subs, behPath, eegPath);
end
end

%% Validation functions
function validate_glm_lvl(x)
if ~ismember(x, [1 2])
    errorMsg = 'Invalid GLM level specified; must be 1 or 2.';
    error(errorMsg);
end
end

%% Other helper functions
function pred = get_predictor_RDM(subs, behPath, eegPath)
%{
        1. BEHAVIORAL RATINGS FROM fMRI
        °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
%}
fmriSubIdx = [1:4 6:subs(end)];
fmriSubs = length(fmriSubIdx);

% define params.
runs = 7;
cats = 2;  % good & bad images
imgsPerCat = 50;
goodRatings = zeros(imgsPerCat, fmriSubs*runs);
badRatings = zeros(imgsPerCat, fmriSubs*runs);
colCount = 0;  % initialize column count for rating arrays

for iSub = 1:fmriSubs
    behSubID = ['sub-' sprintf('%02d', fmriSubIdx(iSub))];

    % get behavioral `.mat` files
    allFiles = dir(fullfile(behPath, behSubID, '*.mat'));

    for iRun = 1:runs
        file = load(fullfile(behPath, behSubID, allFiles(iRun).name));
        colCount = colCount + 1;  % advance column count

        % loop over images
        for iCat = 1:cats
            for iImg = 1:imgsPerCat
                % get correct rating (and handle annoying exceptions 🤢)
                if iSub == 7 && iRun == 7 && iCat == 1 && iImg == 13
                    rating = 0;
                elseif iSub == 10 && iRun == 4 && iCat == 2 && iImg == 19
                    rating = 0;
                else
                    rating = file.dat.resp( ...
                        file.dat.randomTrials(:,1) == iCat & ...
                        file.dat.randomTrials(:,2) == iImg ...
                        );
                end

                % linear interpolation to rescale values to 1-to-7 scale
                rating = 1 + (rating - 1) * (7 - 1) / (4 - 1);

                % add rating to array
                if iCat < 2
                    goodRatings(iImg, colCount) = rating;
                else
                    badRatings(iImg, colCount) = rating;
                end
            end
        end
    end
end

% concatenate rating arrays (good images on top)
allRatings = [goodRatings; badRatings];

% row-wise average
avgRatings = mean(allRatings, 2);

% create difference matrix
fmriRatingDiffs = zeros(100, 100);

for i = 1:size(avgRatings, 1)
    for j = 1:size(avgRatings, 1)
        fmriRatingDiffs(i,j) = abs(avgRatings(i) - avgRatings(j));
    end
end

%{
        2. BEHAVIORAL RATINGS FROM EEG
        °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
%}
eegBehRes = load(fullfile(eegPath, 'behavioral_responses.mat'));
eegBehRes = eegBehRes.behavioral_responses;

% calculate absolute pairwise differences between each image's avg. rating
eegRatingDiffs = zeros(size(fmriRatingDiffs));  % preallocate

for i = 1:size(eegBehRes.rating, 2)
    for j = 1:size(eegBehRes.rating, 2)
        eegRatingDiffs(i,j) = abs( ...
            mean(eegBehRes.rating(:,i)) - mean(eegBehRes.rating(:,j)) ...
            );
    end
end

%{
        3. ADD THE RDMs AND AVERAGE THEM
        °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
%}
pred = (fmriRatingDiffs + eegRatingDiffs) ./ 2;
end
