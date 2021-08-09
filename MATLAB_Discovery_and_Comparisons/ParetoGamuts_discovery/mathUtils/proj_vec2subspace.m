function proj_ss_v = proj_vec2subspace(col_v, subspace)
    % project vector v onto the desired subspace
    % assume ambient space is N dimensional, and subspace is n<N dimensional
    % @param v -- N dimensional COLUMN vector to be projected
    % @param subspace -- nxN matrix defining the desired subspace
    % @return proj_ss_v -- Nx1 COLUMN vector projected onto subspace
    
    nullspace = null(subspace);
    nullspace = orth(nullspace); % don't double count any components
    
    % subtract all components of v that are orthogonal to (in null space of)subspace
    % each column of nullspace is an orthogonal direction
    proj_ss_v = col_v;
    for col=1:size(nullspace, 2)
        null_dir = nullspace(:, col);
        proj_ss_v = proj_ss_v - proj_vec2vec(col_v, null_dir);
    end
end