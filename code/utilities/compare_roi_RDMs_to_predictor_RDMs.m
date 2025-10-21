function d = compare_roi_RDMs_to_predictor_RDMs(d, cfg)
%  COMPARE_ROI_RDMS_TO_PREDICTOR_RDMS Brief summary of this function.
%
% Detailed explanation of this function.
% evaluate input
if ~isfield(cfg, 'RDM_to_partial_out'); cfg.RDM_to_partial_out = {'typical_late', 'control_late'}; end
if ~isfield(cfg, 'correlation_type'); cfg.correlation_type = 'pearson';end
if ~isfield(cfg, 'plot_rdm'); cfg.plot_rdm = false;end
if ~isfield(cfg, 'permutation_test'); cfg.permutation_test = false;end
if ~isfield(cfg, 'n_permutations'); cfg.n_permutations = 10000;end
if ~isfield(cfg, 'add_legend'); cfg.add_legend = false;end
if ~isfield(cfg, 'show_single_cate'); cfg.show_single_cate = false;end
if ~isfield(cfg, 'permutation_type'); cfg.permutation_type = 'row_col_shuffle_ref';end
if ~isfield(cfg, 'order_predictors'); cfg.order_predictors = false;end
if ~isfield(cfg, 'partial_cor'); cfg.partial_cor = true;end
if ~isfield(cfg, 'save_name'); cfg.save_name = 'compare_roi_RDMs_to_predictor_RDMs';end
if ~isfield(cfg, 'xaxis_labels'); cfg.xaxis_labels = true;end
cfg.plot_rdm = false;
if ~isfield(cfg, 'plot_type'); cfg.plot_type = 'bar';end
if ~isfield(cfg, 'dnns'); cfg.dnns = {cfg.dnn};end
if ~isfield(cfg, 'scatter_in_violin'); cfg.scatter_in_violin = 0;end
if cfg.scatter_in_violin == 0
    violin_type = 'full';
elseif cfg.scatter_in_violin == 1
    violin_type = 'half';
end
if ~isfield(cfg, 'ylim'); cfg.ylim = [-0.4, 0.4];end
if ~isfield(cfg, 'task_plotting'); cfg.plott_gap = 0;end
if ~isfield(cfg, 'plotting_predictors'); cfg.plotting_predictors = 1:numel(cfg.predictor_RDMs);end
% get variable attributes
colors = zeros(numel(cfg.RDM_to_partial_out),3);
short_names = cell(1, numel(cfg.RDM_to_partial_out));
for var = 1:numel(cfg.RDM_to_partial_out)
    if strcmp(cfg.RDM_to_partial_out{var}, 'typical_late')
        short_names{var} = 'Typcial drawing';
        colors(var,:) = [1, 0, 1];
    elseif strcmp(cfg.RDM_to_partial_out{var}, 'control_late')
        short_names{var} = 'Control drawing';
        colors(var,:) = [.8, .8, .8];
    elseif strcmp(cfg.RDM_to_partial_out{var}, 'photos_late')
        short_names{var} = 'Photos';
        colors(var,:) = [.4, .9, 1];
    end
end
% prepare figure
if cfg.plotting
    figure;
    hold on
    previous_x_pos = 0;
end
% prepare random permutation (each roi and category should have the same
% random samplings)
if cfg.permutation_test
    % generate permutated subjects list
    rng('default')
    rng(0)
    random_seqs = cell(numel(cfg.RDM_to_partial_out), cfg.n_permutations);
    for i = 1:cfg.n_permutations
        if ismember(cfg.permutation_type, {'row_col_shuffle_ref',...
                'row_col_shuffle_pred', 'row_col_shuffle_pred_all',...
                'row_col_shuffle_pred_plus_reoder'})
            % get random sequence of rows and columns
            for ii = 1:numel(cfg.RDM_to_partial_out)
                random_seqs{ii,i} = randperm(cfg.n);
            end
        elseif strcmp(cfg.permutation_type, 'sign_flip_ref')
            % get random rows and columns to flip sign
            random_seqs{1,i} = randperm(cfg.n, randi([1, cfg.n]));
        end
    end
end

% loop through rois
for roi_i = 1:numel(cfg.rois_of_interest)
    roi = char(cfg.rois_of_interest(roi_i));
    cfg.currentroi = roi;
    % loop through categories
    for cate_num = 1:numel(cfg.categories)
        category = char(cfg.categories{cate_num});
        RDMs = [];
        % get roi RDM names
        all_ref_names = {d.([category,'_RDM']).(cfg.ISC_type).name};
        % get roi RDM
        ref_idx = find(strcmp(all_ref_names, roi));
        if ~isempty(ref_idx)
            ref_RDM = d.([category,'_RDM']).(cfg.ISC_type)(ref_idx);
        else
            tem_preditor_RDMs = cfg.predictor_RDMs; % store predictors temporally
            cfg.predictor_RDMs = {roi};
            RDMs = d.([category,'_RDM']).(cfg.ISC_type)(1); % fill with some RDM as placeholder
            labels = {RDMs.name};
            [ref_RDM, cfg.labels] = evaluate_predictor_RDMs(d, RDMs, labels, cfg, category);
            cfg.predictor_RDMs = tem_preditor_RDMs; % write back predictors
            ref_RDM = ref_RDM(2:end); % remove placeholder
            cfg.labels = cfg.labels{2:end};
        end
        % get canditate/predictor RDMs
        RDMs.name = ref_RDM.name;
        RDMs.color = ref_RDM.color;
        RDMs.RDM = ref_RDM.RDM;
        labels = {ref_RDM.name};
        [RDMs, cfg.labels] = evaluate_predictor_RDMs(d, RDMs, labels, cfg, category);
        % make a cell that holds all predictor RDM structs
        for field = 1:numel({RDMs.name})
            RDMs(field).name = char(cfg.labels{field}); % give it a comprehensive name
        end
        % partial correlation
        if cfg.partial_cor
            [~, r_mat, ~, cfg] = partial_cor_RDM(cfg, RDMs);
        else
            [~, r_mat, ~] = cor_RDM(RDMs,cfg);
        end
        % store results in table
        res_table = table;
        res_table.name = cfg.labels(2:end)';
        res_table.r_val = r_mat(2:end, 1);

        % store in data struct
        d.compare_roi_to_predictor.(roi).(category) = res_table;
        % make random permutations
        if cfg.permutation_test
            permutation_RDMs = RDMs;
            perm_r_mat = zeros(height(r_mat)-1, cfg.n_permutations);
            for perm = 1:cfg.n_permutations
                % randomize RDM
                if strcmp(cfg.permutation_type, 'row_col_shuffle_ref')
                    % shuffle the order of rows and columns
                    ref_RDM = RDMs(1);
                    ref_RDM.RDM = ref_RDM.RDM(random_seqs{1,perm}, random_seqs{1,perm});
                    % replace reference RDM by permutated RDM in RDMs struct
                    permutation_RDMs(1) = ref_RDM;
                elseif ismember(cfg.permutation_type, ...
                        {'row_col_shuffle_pred_all','row_col_shuffle_pred_plus_reoder'})
                    % re-order predictors
                    if strcmp(cfg.permutation_type, 'row_col_shuffle_pred_plus_reoder')
                        permutation_RDMs(2:end) = permutation_RDMs(randperm(numel(cfg.RDM_to_partial_out))+1);
                    end
                    % shuffle ros and columns in predictors
                    for pred = 1:numel(cfg.RDM_to_partial_out)
                        % get predictor RDM
                        pred_RDM = permutation_RDMs(pred+1).RDM;
                        % replace predictor RDM by permutated RDM in RDMs struct
                        permutation_RDMs(pred+1).RDM = ...
                            pred_RDM(random_seqs{pred, perm}, random_seqs{pred, perm});
                    end
                elseif strcmp(cfg.permutation_type, 'sign_flip_ref')
                    ref_RDM_rand = ref_RDM;
                    % flip sign of row
                    ref_RDM_rand.RDM(random_seqs{1,perm},:) = -(ref_RDM.RDM(random_seqs{1,perm},:));
                    % flip sign of columns
                    ref_RDM_rand.RDM(:,random_seqs{1,perm}) = -(ref_RDM.RDM(:,random_seqs{1,perm}));
                    % replace reference RDM by permutated RDM in RDMs struct
                    permutation_RDMs(1) = ref_RDM_rand;
                end
                % run partial correlation
                if ~strcmp(cfg.permutation_type, 'row_col_shuffle_pred')
                    [~, r_mat, ~, ~] = partial_cor_RDM(cfg, permutation_RDMs);
                    perm_r_mat(1:end, perm) = r_mat(2:end, 1);
                else
                    % if permutation of only one predictor while leaving the one to partial out
                    % intact we have to loop through the predictors and do correaltions
                    % seperately
                    % shuffle ros and columns in predictors
                    for pred = 1:numel(cfg.RDM_to_partial_out)
                        % get predictor RDM
                        permutation_RDMs = RDMs;
                        pred_RDM = permutation_RDMs(pred+1).RDM;
                        % replace predictor RDM by permutated RDM in RDMs struct
                        permutation_RDMs(pred+1).RDM = pred_RDM(random_seqs{pred, perm}, random_seqs{pred, perm});
                        [~, r_mat, ~, ~] = partial_cor_RDM(cfg, permutation_RDMs);
                        perm_r_mat(pred, perm) = r_mat(pred+1, 1);
                    end
                end
                if mod(perm/cfg.n_permutations, 0.1) == 0
                    disp([num2str((perm/cfg.n_permutations)*100),...
                        '% of permutations of ', roi, ' ', category, ' is done'])
                end
            end
            d.compare_roi_to_predictor.permutation_test.(roi).(category) = perm_r_mat;
        end

        % runtime control
        disp(['Compare inter-subject RDM of ', roi, ' with partial correaltion of predictor RDMs - ', category])

        % store in data struct
        d.compare_task_to_predictor.(roi).(category) = res_table;
    end
    % average categories
    res_table.r_val_cate1 = d.compare_roi_to_predictor.(roi).(cfg.categories{1}).r_val;
    res_table.r_val_cate2 = d.compare_roi_to_predictor.(roi).(cfg.categories{2}).r_val;
    res_table.r_val = (d.compare_roi_to_predictor.(roi).(cfg.categories{1}).r_val +...
        d.compare_roi_to_predictor.(roi).(cfg.categories{2}).r_val)/2;
    d.compare_roi_to_predictor.(roi).category_average = res_table;
    % get confidence intervals
    if cfg.permutation_test
        % get p values of random permutation
        perm_r_mat = (d.compare_roi_to_predictor.permutation_test.(roi).(cfg.categories{1}) +...
            d.compare_roi_to_predictor.permutation_test.(roi).(cfg.categories{2}))/2;
        d.compare_roi_to_predictor.permutation_test.(roi).category_average = perm_r_mat;
    end
    % plotting
    if cfg.plotting

        % filter predictors to plot
        if numel(cfg.plotting_predictors) ~= numel(cfg.predictor_RDMs)
            res_table = res_table(cfg.plotting_predictors, :);
        end

        if ~isempty(short_names) && numel(cfg.plotting_predictors) ~= numel(short_names)
            short_names = short_names(cfg.plotting_predictors);
            colors = colors(cfg.plotting_predictors, :);
        end

        if numel(cfg.plotting_predictors) ~= height(perm_r_mat)
            perm_r_mat = perm_r_mat(cfg.plotting_predictors, :);
        end

        % loop through res_table and add according variables
        if cfg.plotting_predictors == 1
            if strcmp(roi, 'V1')
                res_table.color_R = .97;
                res_table.color_G = .41;
                res_table.color_B = .92;
            elseif strcmp(roi, 'LOC')
                res_table.color_R = 1;
                res_table.color_G = 0;
                res_table.color_B = 1;
            elseif strcmp(roi, 'PPA')
                res_table.color_R = .73;
                res_table.color_G = .02;
                res_table.color_B = .73;
            elseif strcmp(roi, 'TOS')
                res_table.color_R = .47;
                res_table.color_G = .03;
                res_table.color_B = .45;
            elseif strcmp(roi, 'LPFC')
                res_table.color_R = .24;
                res_table.color_G = .04;
                res_table.color_B = .3;
            else
                res_table.color_R = 1;
                res_table.color_G = 0;
                res_table.color_B = 1;
            end

        else
            clr = cool(numel(RDMs));
            for row = 1:height(res_table)
                % colors
                if ~isempty(colors)
                    res_table.color_R(row) = colors(row,1);
                    res_table.color_G(row) = colors(row,2);
                    res_table.color_B(row) = colors(row,3);
                else % if no colors specified make the bars grey
                    res_table.color_R(row) = clr(row,1);
                    res_table.color_G(row) = clr(row,2);
                    res_table.color_B(row) = clr(row,3);
                end
            end
        end
        for row = 1:height(res_table)
            % short names
            if ~isempty(short_names)
                res_table.short_names{row} = short_names{row};
            else
                res_table.short_names{row} = cfg.labels{row+1};
            end

            % p val
            if cfg.permutation_test
                % get p value from random permutation (one-sided test aginst
                % permutation distribution)
                p_value = sum(perm_r_mat(row,:) >= res_table.r_val(row)) / cfg.n_permutations;
                res_table.p_val(row) = p_value;
                % get p values for categories
                res_table.p_val_cate1(row) = ...
                    sum(d.compare_roi_to_predictor.permutation_test.(roi).(cfg.categories{1})(row,:)...
                    >= res_table.r_val_cate1(row)) / cfg.n_permutations;
                res_table.p_val_cate2(row) = ...
                    sum(d.compare_roi_to_predictor.permutation_test.(roi).(cfg.categories{2})(row,:)...
                    >= res_table.r_val_cate2(row)) / cfg.n_permutations;
                res_table.ci_upper(row) = res_table.r_val(row) - prctile(perm_r_mat(row,:), 5);
                res_table.ci_lower(row) = res_table.r_val(row) - prctile(perm_r_mat(row,:), 95);
            else
                % get p values from r values
                N = nchoosek(cfg.n, 2);
                res_table.p_val(row) = r2p(res_table.r_val(row), N);
                res_table.p_val_cate1(row) = r2p(res_table.r_val_cate1(row), N);
                res_table.p_val_cate2(row) = r2p(res_table.r_val_cate2(row), N);
            end
        end

        % store in data struct
        d.compare_task_to_predictor.(roi).category_average = res_table;
        % order row based on r values
        if cfg.order_predictors
            res_table = sortrows(res_table, 'r_val', 'descend');
        end

        % make plot
        if strcmp(cfg.plot_type, 'violin')
            for xiPos = 1:height(res_table)
                current_x_pos = previous_x_pos + xiPos;
                if cfg.permutation_test
                    Y = {perm_r_mat(xiPos,:)'};
                end

                % make violin plot
                currentColor = [res_table.color_R(xiPos), res_table.color_G(xiPos), res_table.color_B(xiPos)];
                mainHandles(current_x_pos) = daviolinplot(Y,...
                    'color', currentColor,...
                    'violin', violin_type, 'violinalpha', 0.2,...
                    'scatter',cfg.scatter_in_violin,'scatteralpha',0.2,'jitter',1,'scattercolors', 'same', 'scattersize', 5,...
                    'box', 0,...
                    'outliers',0);
                mainHandles(current_x_pos).ds.Vertices(:, 1) = mainHandles(current_x_pos).ds.Vertices(:, 1) + current_x_pos - 1;
                mainHandles(current_x_pos).ds.EdgeColor = currentColor;
                mainHandles(current_x_pos).ds.EdgeAlpha = 0.5;
                mainHandles(current_x_pos).ds.LineWidth = 2;
                if cfg.scatter_in_violin == 1
                    mainHandles(current_x_pos).sc.XData = mainHandles(current_x_pos).sc.XData + current_x_pos - 1.05;
                    mainHandles(current_x_pos).sc.MarkerEdgeColor = currentColor;
                    mainHandles(current_x_pos).sc.MarkerEdgeAlpha = 0.2;
                end
            end
        elseif strcmp(cfg.plot_type, 'bar')
            for xiPos = 1:height(res_table)
                current_x_pos = previous_x_pos + xiPos;
                % Draw individual bar
                barColor = [res_table.color_R(xiPos),res_table.color_G(xiPos),res_table.color_B(xiPos)];
                mainHandles(current_x_pos) = bar(current_x_pos, res_table.r_val(xiPos),...
                    'FaceColor', barColor, 'EdgeColor', 'k');
            end
        end
        for xiPos = 1:height(res_table)
            current_x_pos = previous_x_pos + xiPos;
            if strcmp(cfg.plot_type, 'bar')
                if contains(res_table.name(xiPos), 'Originhal')
                    % Apply hatch only to this bar
                    hatchfill2(mainHandles(current_x_pos), 'HatchAngle', 45, ...
                        'HatchColor', 'k', ...
                        'HatchLineWidth', 1);
                end
            end

            if cfg.permutation_test
                if strcmp(cfg.plot_type, 'bar')
                    % add confidence interval if available
                    if ismember('ci_lower', res_table.Properties.VariableNames)
                        r_val = res_table.r_val(xiPos);
                        errorHandles(current_x_pos) = errorbar(current_x_pos,...
                            r_val, r_val-res_table.ci_lower(xiPos),...
                            res_table.ci_upper(xiPos)-r_val, 'k', 'LineWidth', 1);  % Error bars
                    end
                elseif strcmp(cfg.plot_type, 'violin')
                    % plot oberseved mean r
                    plot([current_x_pos-0.15,current_x_pos+0.15], [res_table.r_val(xiPos), res_table.r_val(xiPos)],...
                        'Color', [res_table.color_R(xiPos), res_table.color_G(xiPos), res_table.color_B(xiPos)], 'LineWidth',3);
                end
            end
            % add marks for single category
            if cfg.show_single_cate
                if cfg.exp_num == 1
                    cate_mark1 = 'B';
                    cate_mark2 = 'K';
                elseif cfg.exp_num == 2
                    cate_mark1 = 'B';
                    cate_mark2 = 'L';
                end
                text(current_x_pos-0.2, res_table.r_val_cate1(xiPos), cate_mark1,...
                    'HorizontalAlignment', 'center', 'FontSize', 5, 'FontWeight', 'bold');
                text(current_x_pos-0.2, res_table.r_val_cate2(xiPos), cate_mark2,...
                    'HorizontalAlignment', 'center', 'FontSize', 5, 'FontWeight', 'bold');
            end
        end
    end
    % make cap between reference RDMs
    previous_x_pos = current_x_pos + cfg.plott_gap;
end

if cfg.plotting

    % collect p values
    all_p_vals = nan(height(res_table), numel(cfg.rois_of_interest));
    for roi_i = 1:numel(cfg.rois_of_interest)
        roi = char(cfg.rois_of_interest(roi_i));
        all_p_vals(:, roi_i) = d.compare_task_to_predictor.(roi).category_average.p_val;
    end

    % get asterisks
    all_asterisks = cell(height(res_table), numel(cfg.rois_of_interest));
    for i_pred = 1:height(res_table)
        % do fdr correction
        [~, ~, ~, fdr_pval] = fdr_bh(all_p_vals(i_pred, :));
        all_asterisks(i_pred, :) = pval2asterisks(fdr_pval, 'none');

        % write back adjusted p values and print them
        disp([newline, newline])
        disp(char(res_table.name(i_pred)))
        for roi_i = 1:numel(cfg.rois_of_interest)
            roi = char(cfg.rois_of_interest(roi_i));
            d.compare_task_to_predictor.(roi).category_average.p_val(i_pred) = fdr_pval(roi_i);
            disp(['FDR corrected p value for ', roi, ': ', num2str(fdr_pval(roi_i)),...
                ' ', char(all_asterisks(i_pred, roi_i))])
            disp(['R value for ', roi, ': ', ...
                num2str(d.compare_task_to_predictor.(roi).category_average.r_val(i_pred))])
        end
    end

    % plot asterisks
    ast_vec = reshape(all_asterisks, 1, []);
    y_pos = zeros(1, length(mainHandles));
    count = 0;
    for i_bar = 1:length(mainHandles)
        if strcmp(cfg.plot_type, 'bar')
            if ~isgraphics(mainHandles(i_bar))
                continue
            end

            % get y position based on error bars
            y_pos(i_bar) = errorHandles(i_bar).YPositiveDelta + errorHandles(i_bar).YData + 0.08;

        elseif strcmp(cfg.plot_type, 'violin')
            if ~isgraphics(mainHandles(i_bar).ds)
                continue
            end
            y_pos(i_bar) = max(mainHandles(i_bar).ds.Vertices(:, 2)) + 0.08;

        end
        count = count + 1;
        astHandle{count} = text(i_bar, y_pos(i_bar), ast_vec{count}, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 20);
    end

    for i_ast = 1:length(mainHandles)
        astHandle{i_ast}.Position(2) = max(y_pos);
    end

    % get aesthetics
    hold off
    if cfg.partial_cor
        ylabel(['Partial correlation [r]', newline]);
    else
        ylabel([cfg.correlation_type, ' correlation [r]', newline]);
    end
    %title('Compare reference RDM with predictors')
    if isfield(cfg, 'plot_type')
        ylim(cfg.ylim)
    else
        ylim([-0.1, max(res_table.r_val) + 0.1])
    end
    xlim([0, previous_x_pos + 1])
    set(gca, 'LineWidth', 2, 'FontName', cfg.FontName, 'FontSize', cfg.FontSize, 'FontWeight', 'bold')
    ax = gca;
    ax.Box = 'off';
    yline(0, 'LineWidth', 2, 'Color', 'k');
    % get labels
    if cfg.xaxis_labels
        xticks(ceil(length(cfg.plotting_predictors)/2):...
            length(cfg.plotting_predictors)+cfg.plott_gap:...
            (length(cfg.plotting_predictors)+2)*length(cfg.rois_of_interest));
        xticklabels(cfg.rois_of_interest);
        xtickangle(45);
    else
        xticklabels([]);
        ax.XColor = 'none';
    end
    % add legend to last plot
    if cfg.add_legend
        legend(res_table.short_names, 'Location','northeastoutside');
    end
    % saving
    fig_path = fullfile(pwd, 'figures', ['exp_', num2str(cfg.exp_num)], 'compare_roi_RDMs_to_predictor_RDMs');
    save_plot(cfg.save_name, fig_path)
end
end