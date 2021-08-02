function  averageNeighbors(buffer,idx)
    i = idx(1);
    j = idx(2);
    p = buffer.Buff(i,j);
    neighboringIndices = [i+1 j-1; i+1 j; i+1 j+1; i j-1; i j+1; i-1 j-1; i-1 j;i-1 j+1];
    meanD = zeros(1,buffer.rD);
    meanF = zeros(1,buffer.rd);
    counter = 0;
    for k = 1:size(neighboringIndices,1)
        idx = neighboringIndices(k,:);
        neighborBuff = buffer.Buff(idx(1),idx(2));

        if buffer.emptyCells(idx(1),idx(2)) == 0
                            meanD = meanD + neighborBuff.minD;
                            meanF = meanF + neighborBuff.minF;
                            counter = counter + 1;
                            selectedPatch = buffer.selectedPatches(idx(1),idx(2));
        end
    end
                        buffer.Buff(i,j).minD = meanD/counter;
                        buffer.Buff(i,j).minF = meanF/counter;
                        buffer.selectedPatches(i,j) = selectedPatch;
end

