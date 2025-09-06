function perm_r_vals = doPermutations(cfg, vec_RDM1, vec_RDM2, vec_RDM3)

% preallocate
ref_RDM = squareform(vec_RDM1);
perm_r_vals = zeros(1, cfg.n_permutations);
rng('default')
rng(1)

for perm = 1:cfg.n_permutations

    % randomize RDM
    % shuffle the order of rows and columns
    shuffledRDM = ref_RDM(cfg.random_seqs{perm}, cfg.random_seqs{perm});

    % partial correlation
    perm_r_vals(perm) = partialcorr(squareform(shuffledRDM)', vec_RDM2', vec_RDM3',...
        'Tail', 'right', 'Type', cfg.partial_correlation_type);

end
end