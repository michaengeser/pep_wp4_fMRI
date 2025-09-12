% --- User inputs ---
vol_fn = 'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\GitHub\pep_wp4_fMRI\ISCtoolbox\wholeBrainAverageTresholded_p001_cluster30_edge.nii';

% --- Load volumes ---
data = load_untouch_nii(vol_fn);
dataImg = data.img;
mask = dataImg > 0; % convert to logical

% --- Label clusters (26-connectivity) ---
CC = bwconncomp(mask, 18);
nClusters = CC.NumObjects;
fprintf('Found %d clusters\n', nClusters);

results = struct();

for c = 1:nClusters
    voxIdx = CC.PixelIdxList{c};              % linear indices in mask
    [i, j, k] = ind2sub(size(mask), voxIdx);    % voxel subscripts

    % Peak (voxel with maximum statistic within this cluster)
    maxR = max(dataImg(voxIdx));
    linearMax = find(dataImg == maxR);
    [imax, jmax, kmax] = ind2sub(size(mask), linearMax);

    % Center of mass (voxel-wise average; can weight by stat if desired)
    mean_sub = round(mean([i j k], 1));

    results(c).voxIdx = voxIdx;
    results(c).peak_mni =  [imax, jmax, kmax];
    results(c).centroid_mni = mean_sub;
end

% --- User inputs: atlas NIfTI and label list (cell array 'labels') ---
atlas_fn = fullfile(pwd, '..', 'MNI_ROIs', 'waal.nii');      % atlas aligned to same space (MNI)
label_file = fullfile(pwd, '..', 'MNI_ROIs', 'aal_labels.txt'); % optional csv with label numbers & names

Vatla = spm_vol(atlas_fn);
atlas_img = spm_read_vols(Vatla);

% If you have a CSV file mapping values->names, load it
% Expect CSV with columns: index,label
labels = {};
if exist(label_file,'file')
    T = readtable(label_file,'ReadVariableNames',false);
    % e.g., T.Var1 = index, T.Var2 = label name
    for r=1:height(T)
        labels{T.Var1(r)} = T.Var2{r};
    end
end

% Loop over clusters and get atlas value at peak and centroid
for c = 1:nClusters
    % get coords
    ix = results(c).peak_mni(1); iy = results(c).peak_mni(2); iz = results(c).peak_mni(3);

    % Check bounds
    if ix>=1 && ix<=size(atlas_img,1) && iy>=1 && iy<=size(atlas_img,2)...
            && iz>=1 && iz<=size(atlas_img,3)
        atlas_val_peak = atlas_img(ix,iy,iz);
    else
        atlas_val_peak = 0;
    end

    results(c).atlas_val_peak = atlas_val_peak;
    if atlas_val_peak>0 && atlas_val_peak<=numel(labels) && ~isempty(labels{atlas_val_peak})
        results(c).atlas_label_peak = labels{atlas_val_peak};
    else
        results(c).atlas_label_peak = sprintf('Label_%g', atlas_val_peak);
    end

    % same for centroid
        % get coords
    ix_cent = results(c).centroid_mni(1); iy_cent = results(c).centroid_mni(2);...
        iz_cent = results(c).centroid_mni(3);
    if ix_cent>=1 && ix_cent<=size(atlas_img,1) && iy_cent>=1 && iy_cent<=size(atlas_img,2)...
            && iz_cent>=1 && iz_cent<=size(atlas_img,3)
        val_cent = atlas_img(ix_cent, iy_cent, iz_cent);
    else
        val_cent = 0;
    end
    
    results(c).atlas_val_centroid = val_cent;
    if val_cent>0 && val_cent<=numel(labels) && ~isempty(labels{val_cent})
        results(c).atlas_label_centroid = labels{val_cent};
    else
        results(c).atlas_label_centroid = sprintf('Label_%g', val_cent);
    end

    fprintf('Cluster %d: peak label = %s (atlas value %g)\n', ...
        c, results(c).atlas_label_peak, results(c).atlas_val_peak);
end
