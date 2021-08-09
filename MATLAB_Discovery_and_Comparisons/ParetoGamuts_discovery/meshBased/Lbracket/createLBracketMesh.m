function m = createLBracketMesh(length, width, thickness, angle)
    %BRACKET create the mesh geometry
    %  for now, assume both arms of L bracket are symmetric
    %  
    %    v1 -- v2
    %     |     |
    %     |     |
    %     |    v4 ------ v6
    %     |              |
    %     v3 ---------- v5    same order for far side, verts v7-12
    %
    %   length = length of arms, eg from v3-v1 or v5-v3
    %   width = distance between front and back (eg between v3-v9)
    %   thickness = distance between eg v1 and v2
    %   angle = angle between outer arm plates; given in radians 
    
%     addpath('../../mathUtils/');
    
    ref_pt = [0,0,0]; % reference point for all other measurements
    
    %vertices = computeVertices_90(ref_pt, width, length, thickness);
    vertices = computeVertices_angle(ref_pt, width, length, thickness, angle);
    
    faces = [
        1, 4, 2;    % near side
        4, 1, 3;
        4, 3, 5;
        4, 5, 6;
        7, 8, 10;   % far side
        7, 10, 9;
        9, 10, 11;
        10, 12, 11;
        1, 7, 3;    % back
        3, 7, 9;
        10, 8, 2;   % front (inner bracket face)
        10, 2, 4;
        11, 12, 6;  % front (edge)
        11, 6, 5;
        9, 11, 5;   % bottom
        9, 5, 3;    
        8, 7, 1;    % top (edge)
        8, 1, 2;
        12, 10, 4;  % top (inner bracket face)
        12, 4, 6;
        ];
    
    m = Mesh(vertices, faces);
end

function vertices = computeVertices_90(ref_pt, width, length, thickness)
    v9 = ref_pt;                        % far lower left
    v3 = ref_pt + [width, 0, 0];        % near lower left
    
    v1 = v3 + [0, 0, length];           % near upper left
    v7 = v9 + [0, 0, length];           % far upper left 
    
    v5 = v3 + [0, length, 0];           % near lower right
    v11 = v9 + [0, length, 0];          % far lower right
    
    v6 = v5 + [0, 0, thickness];        % near center right
    v12 = v11 + [0, 0, thickness];      % far center right
    
    v2 = v1 + [0, thickness, 0];        % near upper center
    v8 = v7 + [0, thickness, 0];        % far upper center
    
    v4 = v3 + [0, thickness, thickness]; % near center center
    v10 = v9 + [0, thickness, thickness]; % far center center
    
    vertices = [v1; v2; v3; v4; v5; v6; v7; v8; v9; v10; v11; v12];
end

function vertices = computeVertices_angle(ref_pt, width, ...
                        length, thickness, angle)
    rotMat = rot([1, 0, 0], angle); % rotation matrix about x axis
    
    offset = [width, 0, 0];
    % define invariant side
    v3 = ref_pt + offset;
    v5 = v3 + [0, length, 0];
    v6 = v5 + [0, 0, thickness];
    
    % define pts to be rotated
    v1 = v5;
    v2 = v1 - [0, 0, thickness];
    
    % rotate the "vertical" beam as desired
    v1 = (rotMat * v1')';
    v2 = (rotMat * v2')';
    
    % compute location of v4; intersection of inner bracket faces
    % we have 1 known pt on each, and direction parallel to outer face
    % only in-plane (yz) so pass last 2 coords
    m_rot = slope2D(v1(2:3), v3(2:3));
    m_horz = slope2D(v5(2:3), v3(2:3)); % should always be 0
    [y,z] = intersectLines2D(v2(2:3), m_rot, v6(2:3), m_horz);

    v4 = [width, y, z];
    
    % extrude into far plane
    v7 = v1 - offset;
    v8 = v2 - offset;
    v9 = v3 - offset;
    v10 = v4 - offset;
    v11 = v5 - offset;
    v12 = v6 - offset;
    
    vertices = [v1; v2; v3; v4; v5; v6; v7; v8; v9; v10; v11; v12];
end
