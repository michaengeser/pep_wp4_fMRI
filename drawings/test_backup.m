% List of open inputs
% Normalise: Write: Images to Write - cfg_files
nrun = X; % enter the number of runs here
jobfile = {'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\GitHub\pep_wp4_fMRI\drawings\test_backup_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(1, nrun);
for crun = 1:nrun
    inputs{1, crun} = MATLAB_CODE_TO_FILL_INPUT; % Normalise: Write: Images to Write - cfg_files
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
