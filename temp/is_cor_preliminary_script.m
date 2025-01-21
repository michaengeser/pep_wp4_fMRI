

fields = fieldnames(RDMstruct);

rdm = RDMstruct.(fields{1});
rdm(eye(size(rdm)) == 1) = 0;
RDM_vec = squareform(rdm);
RDM_mat = zeros(numel(RDM_vec), numel(fields));

for i = 1:numel(fields)

    rdm = RDMstruct.(fields{i});
    rdm(eye(size(rdm)) == 1) = 0;
    RDM_vec = squareform(rdm);
    RDM_mat(:, i) = RDM_vec;

end

meam_rdm = mean(RDM_mat, 2);
meam_rdm = squareform(meam_rdm);
figure; imagesc(meam_rdm, [-1, 1]);

is_rdm = corr(RDM_mat);
figure; imagesc(is_rdm, [-1, 1]);

