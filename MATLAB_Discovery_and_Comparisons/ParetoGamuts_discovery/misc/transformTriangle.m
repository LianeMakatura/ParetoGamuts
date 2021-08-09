rescaledV = V * buffer.numCells + 1;
V_triangle = zeros(size(V));
transformedVertices = zeros(size(V));
colors = zeros(size(V,1),3);
hsv = zeros(size(V,1),3);
clear 'patch'
PIND = buffer.getParetoInd();

for idx = 1:size(V,1)
    vert = rescaledV(idx,:);
    vert = round(vert);

    i = vert(1);
    j = vert(2);
    neighborIndices = [i+1 j-1; i+1 j; i+1 j+1; i j-1; i j+1; i-1 j-1; i-1 j;i-1 j+1];

    if buffer.emptyCells(i, j) == 1
        saturation = -1;
    else
        totalDist = 0;
        counter = 0;
        for k = 1:size(neighborIndices,1) % For each of the neighbors
            r = neighborIndices(k,1);
            c = neighborIndices(k,2);
            if (r < 1 || c < 1 || r > numCells || c > numCells)
                continue
            end    

            if buffer.emptyCells(r,c) == 1
                continue
            end
            
            if max(buffer.Buff(r,c).minD) == 0
                continue;
            end

            diff = norm(buffer.Buff(r,c).minD - buffer.Buff(i,j).minD);
            totalDist = totalDist + diff;
            counter = counter + 1;
        end
        if counter > 0
            saturation = totalDist/counter;
        else
            saturation = -1;
        end
    end
       
    patchColor = buffer.selectedPatches(vert(1),vert(2));
    
    alpha = buffer.getAlpha(vert);
    V_triangle(idx,:) = alpha;
    transformedVertices(idx,:) = alpha;

%     patchColor = patchColor/60;
    colors(idx,:) = [patchColor, patchColor, patchColor];
    value = 1;
    if patchColor == 0
        hue = 0;
    else
        hue = patchColor;
    end
    hsvCol = [hue, saturation, value];
    hsv(idx,:) = hsvCol;
    
end

saturations = hsv(:,2);
saturationsWithoutEmpty = saturations((saturations ~= -1));
minSaturation = min(saturationsWithoutEmpty);
hsv(:,2) = hsv(:,2) - minSaturation;

[Bins,Edges] = histcounts(hsv(:,2));
sortedSats = sort(hsv(:,2));
maxSat = sortedSats(round(1*size(hsv,1)));

for b = 1:size(hsv,1)
    hsv(b,2) = min(hsv(b,2),maxSat);
end
hsv(:,2) = hsv(:,2)/maxSat;
hsv(:,2) = (1-hsv(:,2)).^2;

hueList = unique(hsv(:,1));
deltaHue = 1/size(hueList,1);
for b = 1:size(hsv,1)
    if hsv(b,2) > 1
        hsv(b,2) = 0;
    end
    hu = hsv(b,1);
    
    hsv(b,1) = deltaHue*find(hueList == hsv(b,1));
end

hsv(:,3) = hsv(:,2);

% post process to make everything empty or invalid white
for u = 1:size(hsv,1)
    vert = rescaledV(u,:);
    vert = round(vert);

    i = vert(1);
    j = vert(2);
    
    PIDX = i + (j-1)*buffer.numCells; % Get the flattened index into the buffer of the current point

    if buffer.emptyCells(i,j) == 1
        
        hsv(u,3)= 1;
    elseif ~ismember(PIDX,PIND)
        
        hsv(u,2)= 0;
        hsv(u,3)= 1;
    end
end

converted = hsv2rgb(hsv);
T = T + 1;
figure;
hold on;
plot([0 1 1 0],[0 0 1 1]);
patch('Faces',T,'Vertices',transformedVertices,'FaceVertexCData',converted,'FaceColor','interp', 'LineStyle','none');
T = T-1;
axis off
