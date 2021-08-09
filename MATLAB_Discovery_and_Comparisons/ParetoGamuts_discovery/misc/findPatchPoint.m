function [b] = findPatchPoint(buffer, pIdx, thresh)
    vIdx = buffer.validBufferIndices;
    counter = 0;
    for t = 1:size(vIdx,1)
        idx = vIdx(t);
        b = buffer.Buff(idx);
        if b.bestPatch == pIdx
            counter = counter + 1;
        else
            counter = 0;
        end
        
        if counter == thresh
            break
        end
    end
end

