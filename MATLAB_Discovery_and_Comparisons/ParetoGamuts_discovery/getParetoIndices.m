% takes in an Nxd set of points in performance space, returns logical array
% for each index, that's 1 if non-dominated (should keep it, on Front) or 0
% if dominated (don't keep it)
function result = getParetoIndices(A)

result = ones(size(A,1), 1);
for i = 1:size(A,1)

    P = bsxfun(@minus,A, A(i,:));

    V = zeros(size(P));
    V((P < 0)) = 1;
    if (any (prod(V,2)))
        result(i) = 0;
    end
end
