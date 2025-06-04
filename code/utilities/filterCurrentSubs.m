function d = filterCurrentSubs(d, cfg)

% loop through fields
for iCate = cfg.categories
    iCate = char(iCate);
    for iAna = cfg.analysis_names
        iAna = char(iAna);

        % check which subjects are part of current sample
        dataSubs = d.DNN.(cfg.dnn).(iAna).(iCate).subjects{1, 1};
        currentSubs = nan(1, length(dataSubs));
        for iSub = 1:length(dataSubs)
            currentSubs(iSub) = ismember(dataSubs(iSub), cfg.subNums);
        end
        currentSubs = logical(currentSubs);

        % get RDM and filter for current subjects
        RDM = d.DNN.(cfg.dnn).(iAna).(iCate).subject_mean.RDM;
        d.DNN.(cfg.dnn).(iAna).(iCate).subject_mean.RDM = RDM(currentSubs, currentSubs);
        
    end
end

end