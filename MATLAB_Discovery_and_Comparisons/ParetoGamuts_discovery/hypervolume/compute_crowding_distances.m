function [crowding_distances] = compute_crowding_distances(F)
% [crowding_distances] = compute_crowding_distances(F)
%
% Computes crowding distances given a set of l objective function 
% value vectors. The implementation follows the description of:
%
% 'K. Deb, A. Pratap, S. Argawal, and T. Meyarivan. A Fast and
%  Elitist Multi-Objective Genetic Algorithm: NSGA-II. KanGAL
%  Report No. 2000001, 2001.'
%
% IMPORTANT:
% This function assumes that all solutions of F are non-dominated!
%
% Input:
% - F					- A matrix of M x l, where M is the number
%						  of objectives, and l is the number of
%						  objective function value vectors of the
%						  solutions.
%
% Output:
% - crowding_distances	- a vector of 1 x l with the crowding distances.  
%
% Author: Johannes W. Kruisselbrink
% Last modified: March 17, 2011

	[M, l] = size(F);
	crowding_distances = zeros(1, l);
	for j = 1 : M
		[sort_F, sort_index] = sort(F(j,:));
		crowding_distances(sort_index(1)) = inf;
		crowding_distances(sort_index(end)) = inf;
		for k = 2 : l - 1
			crowding_distances(sort_index(k)) = crowding_distances(sort_index(k)) + (F(j, sort_index(k+1)) - F(j, sort_index(k-1)));
		end
	end

end
