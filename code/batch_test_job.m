%-----------------------------------------------------------------------
% Job saved on 29-Oct-2024 13:42:20 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
%%
matlabbatch{1}.spm.spatial.realign.estwrite.data = {
                                                    {
                                                    'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\GitHub\pep_wp4_fMRI\sourcedata\sub-001\func\sub-001xxxx_task-localizer_bold.nii,1'
                                                    'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\GitHub\pep_wp4_fMRI\sourcedata\sub-001\func\sub-001xxxx_task-rest_bold.nii,1'
                                                    'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\GitHub\pep_wp4_fMRI\sourcedata\sub-001\func\sub-001xxxx_task-scenes_run-10_bold.nii,1'
                                                    'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\GitHub\pep_wp4_fMRI\sourcedata\sub-001\func\sub-001xxxx_task-scenes_run-1_bold.nii,1'
                                                    'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\GitHub\pep_wp4_fMRI\sourcedata\sub-001\func\sub-001xxxx_task-scenes_run-2_bold.nii,1'
                                                    'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\GitHub\pep_wp4_fMRI\sourcedata\sub-001\func\sub-001xxxx_task-scenes_run-3_bold.nii,1'
                                                    'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\GitHub\pep_wp4_fMRI\sourcedata\sub-001\func\sub-001xxxx_task-scenes_run-4_bold.nii,1'
                                                    'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\GitHub\pep_wp4_fMRI\sourcedata\sub-001\func\sub-001xxxx_task-scenes_run-5_bold.nii,1'
                                                    'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\GitHub\pep_wp4_fMRI\sourcedata\sub-001\func\sub-001xxxx_task-scenes_run-6_bold.nii,1'
                                                    'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\GitHub\pep_wp4_fMRI\sourcedata\sub-001\func\sub-001xxxx_task-scenes_run-7_bold.nii,1'
                                                    'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\GitHub\pep_wp4_fMRI\sourcedata\sub-001\func\sub-001xxxx_task-scenes_run-8_bold.nii,1'
                                                    'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\GitHub\pep_wp4_fMRI\sourcedata\sub-001\func\sub-001xxxx_task-scenes_run-9_bold.nii,1'
                                                    }
                                                    }';
%%
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

matlabbatch{2}.spm.spatial.coreg.estimate.ref(1) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
matlabbatch{2}.spm.spatial.coreg.estimate.source = {'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\GitHub\pep_wp4_fMRI\sourcedata\sub-001\anat\sub-001xxxx_T1w.nii,1'};
matlabbatch{2}.spm.spatial.coreg.estimate.other = {''};
matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

matlabbatch{3}.spm.spatial.preproc.channel.vols(1) = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
matlabbatch{3}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{3}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{3}.spm.spatial.preproc.channel.write = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(1).tpm = {'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\MATLAB\spm12\tpm\TPM.nii,1'};
matlabbatch{3}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{3}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{3}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(2).tpm = {'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\MATLAB\spm12\tpm\TPM.nii,2'};
matlabbatch{3}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{3}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{3}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(3).tpm = {'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\MATLAB\spm12\tpm\TPM.nii,3'};
matlabbatch{3}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{3}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{3}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(4).tpm = {'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\MATLAB\spm12\tpm\TPM.nii,4'};
matlabbatch{3}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{3}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{3}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(5).tpm = {'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\MATLAB\spm12\tpm\TPM.nii,5'};
matlabbatch{3}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{3}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{3}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(6).tpm = {'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\MATLAB\spm12\tpm\TPM.nii,6'};
matlabbatch{3}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{3}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{3}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{3}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{3}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{3}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{3}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{3}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{3}.spm.spatial.preproc.warp.write = [0 0];
matlabbatch{3}.spm.spatial.preproc.warp.vox = NaN;
matlabbatch{3}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
                                              NaN NaN NaN];
matlabbatch{4}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: c1 Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
%%
matlabbatch{4}.spm.spatial.normalise.write.subj.resample = {
                                                            'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\GitHub\pep_wp4_fMRI\sourcedata\sub-001\func\sub-001xxxx_task-localizer_bold.nii,1'
                                                            'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\GitHub\pep_wp4_fMRI\sourcedata\sub-001\func\sub-001xxxx_task-rest_bold.nii,1'
                                                            'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\GitHub\pep_wp4_fMRI\sourcedata\sub-001\func\sub-001xxxx_task-scenes_run-10_bold.nii,1'
                                                            'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\GitHub\pep_wp4_fMRI\sourcedata\sub-001\func\sub-001xxxx_task-scenes_run-1_bold.nii,1'
                                                            'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\GitHub\pep_wp4_fMRI\sourcedata\sub-001\func\sub-001xxxx_task-scenes_run-2_bold.nii,1'
                                                            'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\GitHub\pep_wp4_fMRI\sourcedata\sub-001\func\sub-001xxxx_task-scenes_run-3_bold.nii,1'
                                                            'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\GitHub\pep_wp4_fMRI\sourcedata\sub-001\func\sub-001xxxx_task-scenes_run-4_bold.nii,1'
                                                            'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\GitHub\pep_wp4_fMRI\sourcedata\sub-001\func\sub-001xxxx_task-scenes_run-5_bold.nii,1'
                                                            'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\GitHub\pep_wp4_fMRI\sourcedata\sub-001\func\sub-001xxxx_task-scenes_run-6_bold.nii,1'
                                                            'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\GitHub\pep_wp4_fMRI\sourcedata\sub-001\func\sub-001xxxx_task-scenes_run-7_bold.nii,1'
                                                            'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\GitHub\pep_wp4_fMRI\sourcedata\sub-001\func\sub-001xxxx_task-scenes_run-8_bold.nii,1'
                                                            'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\GitHub\pep_wp4_fMRI\sourcedata\sub-001\func\sub-001xxxx_task-scenes_run-9_bold.nii,1'
                                                            };
%%
matlabbatch{4}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                          78 76 85];
matlabbatch{4}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
matlabbatch{4}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{4}.spm.spatial.normalise.write.woptions.prefix = 'w';


matlabbatch{5}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: c1 Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
matlabbatch{5}.spm.spatial.normalise.write.subj.resample = {'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\GitHub\pep_wp4_fMRI\sourcedata\sub-001\anat\sub-001xxxx_T1w.nii,1'};
matlabbatch{5}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                          78 76 85];
matlabbatch{5}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
matlabbatch{5}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{5}.spm.spatial.normalise.write.woptions.prefix = 'w';
