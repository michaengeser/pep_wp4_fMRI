function plot_HCP_parcels_on_inflated(results)

% This function requires the GIFT toolbox from matlab https://de.mathworks.com/matlabcentral/fileexchange/112570-gift
% and Connectome Workbench https://www.humanconnectome.org/software/connectome-workbench

%----------------------------------------------------------
% INPUT FILES
%----------------------------------------------------------
surf_file = fullfile('C:', 'HCP_files', 'S1200.R.inflated_MSMAll.32k_fs_LR.surf.gii'); % Inflated right hemisphere
dlabel_file = fullfile('C:', 'HCP_files', 'Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');

if ~exist(surf_file, 'file') || ~exist(dlabel_file, 'file')
    error(['Inflated brain and label file must be unter C:\HCP_files\ because ', ...
        'Connectome Workbench does not work well with relative paths. ', ...
        'Alternatively, hard code the path to the MNI_ROIs folder on your machine'])
end


%----------------------------------------------------------
% SETTINGS
%----------------------------------------------------------
highlight_parcels = results(8).sig200Parcels; % parcels to highlight
parcelColors = load(fullfile(pwd, 'utilities', 'parcelColors.mat'));
parcelColors = parcelColors.parcelColors;
border_color = [0.2, 0.2, 0.2]; % gray border color

%----------------------------------------------------------
% LOAD SURFACE
%----------------------------------------------------------
surf = gifti(surf_file);
vertices = surf.vertices;
faces = surf.faces;

%----------------------------------------------------------
% EXTRACT RIGHT HEMISPHERE LABELS FROM .dlabel FILE
%----------------------------------------------------------
% Requires wb_command (Connectome Workbench)
disp('Extracting right hemisphere labels from dlabel...');
system(['wb_command -cifti-separate ' dlabel_file ...
    ' COLUMN -label CORTEX_RIGHT temp_R.label.gii']);

label_gii = gifti('temp_R.label.gii');
labels  = double(label_gii.cdata);

%% --- CREATE BASE PLOT ---
disp('Plotting inflated brain...');
figure('Color','w','Position',[100 100 800 600]);
p = patch('Vertices', vertices, 'Faces', faces, ...
    'FaceVertexCData', repmat([0.9 0.9 0.9], size(vertices,1), 1), ...
    'FaceColor', 'interp', 'EdgeColor', 'none');
axis equal off
camlight('headlight');
camlight('right');
lighting gouraud;
material dull;
view([-90,-90])
camlight('left');

hold on;


%% --- DRAW PARCEL BORDERS ---  (This does not realy work...)
% 
% disp('Drawing parcel borders...');
% unique_labels = unique(labels);
% 
% for i = 1:length(unique_labels)
%     this_label = unique_labels(i);
%     mask = (labels == this_label);
% 
%     % find edges where one vertex is inside and one is outside
%     edges = [faces(:,[1 2]); faces(:,[2 3]); faces(:,[3 1])];
%     edges_mask = mask(edges(:,1)) ~= mask(edges(:,2));
%     boundary_edges = edges(edges_mask,:);
% 
%     % remove duplicates
%     boundary_edges = unique(sort(boundary_edges,2),'rows');
% 
%     plot3(vertices(boundary_edges',1), vertices(boundary_edges',2), vertices(boundary_edges',3), ...
%         'Color', border_color, 'LineWidth', 0.8);  % increase line width for clarity
% end


%% --- HIGHLIGHT SELECTED PARCELS ---
disp('Highlighting selected parcels...');
for i = 1:numel(highlight_parcels)
    parcel_id = highlight_parcels(i);
    mask = (labels == parcel_id);
    p = patch('Vertices', vertices, 'Faces', faces(any(mask(faces),2),:), ...
        'FaceColor', parcelColors{i}, 'EdgeColor', 'none', ...
        'FaceAlpha', 0.9);
end
material([0.4 0.3 0.1 0.4]); 

title('Highlighted HCP Parcels on Inflated Right Hemisphere','FontSize',14);


%----------------------------------------------------------
% Helper function: faces2edges
%----------------------------------------------------------
    function edges = faces2edges(faces)
        edges = [faces(:,[1 2]); faces(:,[2 3]); faces(:,[3 1])];
        edges = sort(edges,2);
        edges = unique(edges,'rows');
    end
end