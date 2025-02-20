counter = 0;
for iRun = 1:cfg.nRuns
    for iSub = 1:numel(cfg.subNums)
        subID = sprintf('sub-%0.3d', cfg.subNums(iSub));

        counter = counter +1;
        folderName = fullfile(pwd, 'derivatives', subID, 'timecourses');
        fileName = ['whole_brain_timecourse_kitchen_run_', num2str(iRun),'.mat'];

        %%
        names{counter, 1} = fullfile(folderName, fileName);
        disp(fullfile(folderName, fileName))


    end

end