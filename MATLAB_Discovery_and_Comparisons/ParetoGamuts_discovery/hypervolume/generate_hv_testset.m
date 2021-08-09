function [F_ND] = generate_hv_testset(M, l, testno)
% [F_ND] = generate_hv_testset(M, l, testno)
%
% Generates a test set for hypervolume computation algorithms, being
% a set of l non-dominated objective function vectors (M objectives).
% The testset is both returned as output and writen to a text-file
% "M"D_"l"_points_"testno:.txt.
%
% IMPORTANT:
%   Considers Minimization of the objective function values!
%
% Input:
% - M              - The number of objectives
% - l              - The number of data points
% - testno         - The test number id used within the filename
%                    saving the testset
%
% Output:
% - F_ND           - An M x l matrix of l non-dominated solutions
%
% Author: Johannes W. Kruisselbrink
% Last modified: March 17, 2011

	x = lhsu(M, M, zeros(M,1), ones(M,1))';
	d = 2 * rand(1,M);
	d = 2 * rand() * ones(1,M);

	i = 1;
	F_ND = [];

	location_points = 400;
	while (length(F_ND) < l)
		X = rand(M,i * location_points);
		F = zeros(M,i * location_points);
		for j = 1:i * location_points
			F(:,j) = sum((abs(repmat(X(:,j),1,M) - x)).^repmat(d,M,1));
		end
		non_dominated_front(F)
		F_ND = F(:,non_dominated_front(F))
		i = i + 1;
	end
	F_ND = F_ND(:,1:l);

	filename = [num2str(M), 'D_', num2str(l), '_points_', num2str(testno), '.txt'];
	dlmwrite(filename, F_ND', 'delimiter', '\t');

	if (M == 2)
		clf
		plot(F_ND(1,:), F_ND(2,:), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
		grid on
	elseif (M == 3)
		plot3(F_ND(1,:), F_ND(2,:), F_ND(3,:), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
		grid on
	end

end

