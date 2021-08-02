function adj = getAdjacencyMatrix(T, buffer, V)
    adj = sparse(size(V,1),size(V,1));
    T = T+1;
    for i = 1:size(T,1)
        a = T(i,1);
        p1 = V(a,:);
        b = T(i,2);
        p2 = V(b,:);
        c = T(i,3);
        p3 = V(c,:);
        
        buffIdxA = p1(1) + (p1(2)-1)*buffer.numCells;
        buffIdxB = p2(1) + (p2(2)-1)*buffer.numCells;
        buffIdxC = p3(1) + (p3(2)-1)*buffer.numCells;

        ab = norm(buffer.Buff(buffIdxA).minF - buffer.Buff(buffIdxB).minF);
        bc = norm(buffer.Buff(buffIdxB).minF - buffer.Buff(buffIdxC).minF);
        ac = norm(buffer.Buff(buffIdxA).minF - buffer.Buff(buffIdxC).minF);
        
        
        if ab ~= 0
            
        end
        adj(a,b) = ab;
        adj(b,a) = ab;
        adj(b,c) = bc;
        adj(c,b) = bc;
        adj(a,c) = ac;
        adj(c,a) = ac;
        
    end
%     adj = sparse(adj);
end

