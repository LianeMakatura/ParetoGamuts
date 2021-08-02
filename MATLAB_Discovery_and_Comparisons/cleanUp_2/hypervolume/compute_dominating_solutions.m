function [dominating_solutions, D] = compute_dominating_solutions(F)
% [dominating_solutions, D] = compute_dominating_solutions(F)
%
% Computes for each solution the number of dominating solutions.
%
% Input:
% - F						- A matrix of M x l, where M is the number
%							  of objectives, and l is the number of
%							  objective function value vectors of the
%							  solutions.
%
% Output:
% - dominating_solutions	- a vector of 1 x l which indicates for each
%							  solution the number of dominating solutions.  
% - D						- domination matrix
%
% Author: Johannes W. Kruisselbrink
% Last modified: March 17, 2011

	[M, l] = size(F);
	D = zeros(l,l);
	dominating_solutions = zeros(1,l);
	for i = 1 : l
		for j = 1 : l
			if (dominates(F(:,j), F(:,i)))
				D(j,i) = 1;
				D(i,j) = -1;
				dominating_solutions(1,i) = dominating_solutions(1,i) + 1;
			end
		end
	end

end
