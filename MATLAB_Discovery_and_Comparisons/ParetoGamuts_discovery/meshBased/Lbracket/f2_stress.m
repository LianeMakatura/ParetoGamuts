% Measures the maximum von Mises stress of an LBracket at a given point
% @param x - rD x n list of design points
% @param z - rA x n list of application points
%
%@return stress - 1 x n list of evaluated points

function stress = f2_stress(x,z)
    % handcoded, since no corresponding symbolic function. First necessarily FD test.
    % unpack the block vars into the variable names used by the function
    x1 = x(1, :).';
    x2 = x(2, :).';
    x3 = x(3, :).';
    z1 = z(1, :).';

    % rescale the inputs appropriately
    length = 0.4 + (x1 .* 0.6); % allowed from 0.4-1, 0 nonsense
    width = 0.2 + (x2 .* 0.8); %allowed from 0.2-1, 0 is nonsense
    thickness = 0.1 + (x3.*0.2); %allowed from 0.1-0.3

    % normalize so zSymb between 0-1 gives angle pi/2 to pi/4
    % eps offsets divide by zero errors that happen at pi/4, pi/2
    eps = 1e-3;
    angle = pi/4 + eps + (z1 .* (pi/4 - 2*eps));
    
    % loop over each point and evaluate separately
    numPts = size(x, 2);
    stress = zeros(1, numPts);
    for i=1:numPts
        % measure the maximum stress of the object under a load
        m = createLBracketMesh(length(i, 1), width(i, 1), thickness(i, 1), angle(i, 1));
        fixed_face = 6;
        loaded_face = 3;
        load =  [0;100;0];
        structuralresults = max_stress(m, fixed_face, loaded_face, load);
        s = norm(structuralresults.VonMisesStress, 1000);
        
        % normalize so between 0 and 1
        stress(1, i) = s / 3e4;
    end
end

