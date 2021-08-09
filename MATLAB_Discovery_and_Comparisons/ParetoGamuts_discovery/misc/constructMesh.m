function [ V, T, T_colors, V_metrics, C, colorMapping, retC, vRef ] = constructMesh( buffer, thresh )
    paretoInd = buffer.getParetoInd();

    V = zeros(size(buffer.validBufferIndices,1),2); % V_i will hold the coordinates of the ith vertex in space
    V_metrics = zeros(size(buffer.validBufferIndices,1),buffer.rd+buffer.rD); % V_metrics_i will hold the parameters (design and objective space) of the ith vertex
    V_metrics_everything = zeros(100);
    T = []; % T_i will hold three indices into @V that define a mesh triangle. Bottleneck right now... we don't know how big T will be based off V
    
    C = zeros(size(buffer.Buff,1)^2,1); % C_i will hold the the color value of the ith vertex, big right now for search purposes but will be truncated later
    C_indices = zeros(size(V,1),1); % Indices used to truncate C after flood fill
    
    numCells = buffer.numCells;
    counter = 1;
    for i = 1:numCells %For all rows
        
        range = numCells-i+1; % Only walk as far as the diagonal
        for j = 1:range
            V(counter,:) = [i j]; % Store the indices of the current point
            buffIdx = i + (j-1)*buffer.numCells; % Get the flattened index into the buffer of the current point
            metrics = [buffer.Buff(buffIdx).minD buffer.Buff(buffIdx).minF]; % Get the design/objective parameters for the current point
            C_indices(counter) = buffIdx;
            
%             if ismember(buffIdx,paretoInd)
                V_metrics(counter,:) = metrics; % Write these parameters into @V_metrics 
%             end
            if j == range-1 % We are at the "end" of the row, so there is only one triangle
                tri = [counter counter+1 counter+numCells-i+1]; % Counterclockwise-defined triangle 
                T = [T; tri];
            elseif j < range
                triLeft = [counter counter+1 counter+numCells-i+1];
                triRight = [counter+1 counter+numCells-i+2 counter+numCells-i+1];

                T = [T; triLeft; triRight];
            end
            counter = counter + 1;
        end
    end
    
    vRef = V;
    V = (V-1)/numCells; % Scale all points to be between 0 and 1
    
    % Coloring
    patchID = 1;
    for i = 1:buffer.numCells
        range = buffer.numCells-i+1;
        for j = 1:range
            idx = i + (j-1)*buffer.numCells;
            
            if (sum(buffer.Buff(idx).minF) > 0) % Replace this with proper empty cells check
                if C(idx) == 0% && ismember(idx,paretoInd)
                    C = performColoring(C, buffer, i, j, patchID);
                    patchID = patchID + 1;
                    
                end
            end
        end
    end
    retC = C;
    C = C(C_indices);
    T_colors = zeros(size(T,1),3);
    
    cmap1 = hot(patchID);
    cmap2 = winter(patchID) ;
    % Combine them into one tall colormap.
    combinedColorMap = [cmap1; cmap2];
    % Pick 15 rows at random.
    randomRows = randi(size(combinedColorMap, 1), [6000, 1]);
    % Extract the rows from the combined color map.
    randomColors = combinedColorMap(randomRows, :);
    patchToColor = ceil(255*randomColors);
    patchToColor = randomColors;
%     colorMapping = randomColors;

    for m = 1:size(T,1)
        vertices = T(m,:);
        
        colors = unique(C(vertices));
        verticesDesignSpace = V_metrics(vertices,1:buffer.rD);
        if sum(sum(verticesDesignSpace)) > 0
            d1 = norm(verticesDesignSpace(1,:) - verticesDesignSpace(2,:));
            d2 = norm(verticesDesignSpace(1,:) - verticesDesignSpace(3,:));
            d3 = norm(verticesDesignSpace(2,:) - verticesDesignSpace(3,:));
            if max([d1 d2 d3]) > thresh
                
                T_colors(m,:) = [1 0 0];
            elseif min(buffer.selectedPatches(vertices)) == 0
                T_colors(m,:) = [0 1 0];
            else
               
                color = mean(patchToColor(buffer.selectedPatches(vertices), :));
                T_colors(m,:) = color;
%                 T_colors(m,:) = [0 1 1];
            end
        else
            T_colors(m,:) = [0 0 1];
        end
        
        
        if numel(colors) == 1

        end
    end
    
    % Define the two colormaps.
    cmap1 = hot(patchID);
    cmap2 = winter(patchID) ;
    % Combine them into one tall colormap.
    combinedColorMap = [cmap1; cmap2];
    % Pick 15 rows at random.
    randomRows = randi(size(combinedColorMap, 1), [patchID, 1]);
    % Extract the rows from the combined color map.
    randomColors = combinedColorMap(randomRows, :);
    colorMapping = ceil(255*randomColors);
    colorMapping = randomColors;
%     colorMapping = randomColors;
    colorMapping = [0 0 0; colorMapping];
    T = T - 1; % 0 indexing for C++

end
 
function newC = performColoring(C, buffer, row, col, color)
    numCells = buffer.numCells;
    rootIdx = row + (col-1)*numCells;
    newC = C;
    newC(rootIdx) = color;
    stack = [row col];
    epsilon = .1;
    while size(stack,1) > 0 % While the stack is not empty
        % Assume that all stack members are valid
        subs = stack(1,:); % Get the coordinates of the top of the stack
        i = subs(1);
        j = subs(2);
        index = i + (j-1)*numCells; % Get the buffer index for the top of the stack
        currentD = buffer.Buff(index).minD;

        neighboringIndices = [i+1 j-1; i+1 j; i+1 j+1; i j-1; i j+1; i-1 j-1; i-1 j;i-1 j+1]; % Get the neighbors of the cell
        for k = 1:size(neighboringIndices,1) % For each of the neighbors
            r = neighboringIndices(k,1);
            c = neighboringIndices(k,2);
            if (r < 1 || c < 1 || r > numCells || c > numCells)
                continue
            end
            newIdx = r + (c-1)*numCells; % Get the buffer index of the neighbor
            neighboringD = buffer.Buff(newIdx).minD;
            
            diff = neighboringD-currentD; % Compare the design-space proximity of the original cell and its neighbor
            
            if buffer.emptyCells(r, c) == 0 && newC(newIdx) == 0 && (norm(diff) < epsilon)% || buffer.selectedPatches(newIdx) == buffer.selectedPatches(index))
                
                coords = [r c];
                stack = [stack; coords];
                newC(newIdx) = color;
            end
        end
        
        stack(1,:) = []; % remove top of the stack
    end
end