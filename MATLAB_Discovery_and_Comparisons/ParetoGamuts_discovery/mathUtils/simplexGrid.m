function [S] = simplexGrid(rd, N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    ptsPerDim = ceil(nthroot(N, (rd-1)));
    spacing = (pi/2) / (ptsPerDim-1);
    
%     G = polarGrid(rd, N);
    rng = 0:spacing:pi/2;
    
    % convert to Cartesian coordinates
    % assumes 2d for now
    G = zeros(ptsPerDim, ptsPerDim, rd);
    for i=1:ptsPerDim
        for j=1:ptsPerDim
            phi = [rng(i), rng(j)];
            x = polar2Cartesian(phi, 1);
            G(i,j, :) = x / norm(x, 1);
        end
    end
    surf(G(:, :, 1), G(:, :, 2), G(:, :, 3))
end





