function d = getISCofRespresentation(cfg, d)

% get number of ROI masks
nmasks=numel(cfg.rois);

% loop through categories
for cate_num = 1:numel(cfg.categories)
    category = char(cfg.categories{cate_num});

    % loop through ISC types
    for ISC_type = cfg.ISC_types
        cfg.ISC_type = char(ISC_type);

        % get RDM of ISC type and category
        % loop through ROIs
        for j=1:nmasks
            mask_label=cfg.rois{j};
            mask_label_short = split(mask_label, '.');
            mask_label_short = mask_label_short{1};
            
            % preallocate RDM matrix
            RDMmat = nan(nchoosek(cfg.nTrials/2, 2), cfg.n);

            % loop through subjects
            for iSub = 1:length(cfg.subNums)
                subID = sprintf('sub-%0.3d', cfg.subNums(iSub));
                subID2 = strrep(subID, '-', '');

                % make a matrix with vectorized RDMs
                rdm = squeeze(d.(cfg.ISC_type).(subID2).(mask_label_short).rdm);
                rdm(eye(size(rdm)) == 1) = 0;

                % filter for category
                if strcmp(category, cfg.categories{1})
                    rdm = rdm(1:cfg.nTrials/2, 1:cfg.nTrials/2);
                else
                    rdm = rdm(cfg.nTrials/2 + 1:end, cfg.nTrials/2 + 1:end);
                end
                RDMmat(:, iSub) = squareform(rdm);
            end

            % make IS-RDM
            [~, mat_out, ~] = make_RDM(RDMmat, cfg);
            if cfg.dissimilarity
                median_mat_out = 1 - mat_out;
                median_mat_out(eye(size(median_mat_out)) == 1) = 0;
                medianISC = median(squareform(median_mat_out), 'omitnan');
            else
                median_mat_out = mat_out;
                median_mat_out(eye(size(median_mat_out)) == 1) = 0;
                medianISC = median(squareform(median_mat_out), 'omitnan');
            end

            % store in structure 
            d.([category,'_RDM']).(cfg.ISC_type)(j).RDM = mat_out;
            d.([category,'_RDM']).(cfg.ISC_type)(j).color = [0, 0, 0];
            d.([category,'_RDM']).(cfg.ISC_type)(j).name = (mask_label_short(2:end));
            d.medianISC.(category).(cfg.ISC_type).(mask_label_short(2:end)) = medianISC;
            
        end % mask
    end % isc type
end % category
end