figure;
newC = reshape(RC,100,100);
imagesc(newC);
clear 'patch'
figure;
hold on
T = T + 1;
    for j = 1:length(TC)
        
        color =TC(j,:);
        vertices = T(j,:);
        
        v1 = V(vertices(1),:);
        v2 = V(vertices(2),:);
        v3 = V(vertices(3),:);
        x = [v1(1) v2(1) v3(1)];
        y = [v1(2) v2(2) v3(2)];
        tri = delaunay(x,y);

        
%         triplot(tri,x,y,'Color',color);
 
        patch(x, y, color, 'LineStyle','none');
    end
T = T - 1;
% pause;
f1 = figure;
hold on;
f2 = figure;
hold on;
    for j = 0:size(CM,1)
        
        idx = find(C == j);
    
        if j == 0
            continue
        else
            color = CM(j,:);
            if sum(color) == 3
                color = [.15 .15 .15];
            end
        end
        
        vertices = V_M(idx,:);
        
        figure(f1);
        scatter3(vertices(:,1),vertices(:,2),vertices(:,3), 'MarkerFaceColor', color);
        figure(f2);
        scatter3(vertices(:,rD+1),vertices(:,rD+2),vertices(:,rD+3), 'MarkerFaceColor', color, 'MarkerEdgeColor', color);

    end
% if rD > 3
%     distFlat = pdist(V_M(:,1:rD));
%     disMat = squareform(distFlat);
%     reduced = cmdscale(disMat,3);
%     figure;
%     scatter3(reduced(:,1),reduced(:,2),reduced(:,3));
% end

% T = T - 1 ;
    