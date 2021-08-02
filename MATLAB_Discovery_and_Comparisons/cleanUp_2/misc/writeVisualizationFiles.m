[V, T, TC, V_M, C, CM, RC, vRef] = constructMesh(buffer, .15);

newVM = [];
for i = 1:size(V_M,1)
    line = V_M(i,:);
    if sum(line) ~= 0
        newVM = [newVM; [i line]];
    end
end

newVMDesign = newVM(:,2:rD+1);

onShapeLine = reshape(newVMDesign',1,rD*size(newVM,1));

infoFile = zeros(6,1);
infoFile(1) = buffer.rD;
infoFile(2) = buffer.rd;
infoFile(3) = size(T,1);
infoFile(4) = size(V,1);
infoFile(5) = size(CM,1);
infoFile(6) = size(newVM,1);

filename = '/bike';
if false
dlmwrite(['finalVis' filename '/colors.txt'], CM, ' ');
dlmwrite(['finalVis' filename '/front_v_pos.txt'], V, ' ');
dlmwrite(['finalVis' filename '/front_t.txt'], T, ' ');
dlmwrite(['finalVis' filename '/front_v_label.txt'], TC, ' ');
dlmwrite(['finalVis' filename '/front_v_metrics.txt'], newVM, ' ');
dlmwrite(['finalVis' filename '/info.txt'], infoFile, ' ');
end
% dlmwrite('results/lamp3DVis/onShapeLine.txt', onShapeLine, ',');
% csvwrite('finalVis/onShapeCsv.csv',onShapeLine);
