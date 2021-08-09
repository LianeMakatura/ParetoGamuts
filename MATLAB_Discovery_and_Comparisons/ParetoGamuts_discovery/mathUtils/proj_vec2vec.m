function proj_v_u = proj_vec2vec(u, v)
    % project vector u onto vector v
    proj_v_u = dot(u,v) / (norm(v)^2) * v;
end