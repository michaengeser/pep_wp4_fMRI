function smoothData(cfg)

if ~isfield(cfg, 'nRuns'); cfg.nRuns = 10; end
if ~isfield(cfg, 'smoothKernel'); cfg.smoothKernel = 4; end
if ~isfield(cfg, 'nVols'); cfg.nVols = 152;end


%% Initialize SPM
matlabbatch = [];
spm('defaults', 'fmri');
spm_jobman('initcfg');

%% Define subjects and main path
mainPath = fullfile(pwd, '..');
fmriPath = fullfile(mainPath, 'derivatives');

for iSub = 1:numel(cfg.subNums)
    subID = sprintf('sub-%0.3d', cfg.subNums(iSub));

    disp(['Smoothing subject: ', num2str(cfg.subNums(iSub))])
    
    % for each run, get each scan's `.nii` file
    counter = 0;
    for iRun = 1:cfg.nRuns

        for iVol = 1:cfg.nVols
            counter = counter +1;
            folderName = fullfile(fmriPath, subID, 'func');
            fileName = sprintf('wrsub-%0.3dxxxx_task-scenes_run-%d_bold_%0.5d.nii,1', cfg.subNums(iSub), iRun, iVol);

            %%
            matlabbatch{1}.spm.spatial.smooth.data{counter, 1} = fullfile(folderName, fileName);
        end

    end
    %%
    matlabbatch{1}.spm.spatial.smooth.fwhm = [cfg.smoothKernel, cfg.smoothKernel, cfg.smoothKernel];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';

    % go!
    % repeat 3 times when crashing
    maxRetries = 3;
    for attempt = 1:maxRetries
        try
            spm_jobman('run', matlabbatch);
            break; % Exit loop if successful
        catch ME
            if attempt == maxRetries
                error('Failed after %d attempts: %s', maxRetries, ME.message);
            else
                disp(['Retrying: Attempt ', num2str(attempt)]);
                pause(1); % Small delay before retry
            end
        end
    end

    clear matlabbatch;

end
end