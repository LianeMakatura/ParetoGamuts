function [ref_point] = compute_reference_point(F, offset_fraction)
% [ref_point] = compute_reference_point(F, offset_fraction)
%
% Computes a reference point of the M vectors of l function values
% contained in F. The offset_fraction determines the shift in of the
% reference point as a fraction of the difference between the minimum
% and maximum function value for each direction.
%
% IMPORTANT:
%   Considers Minimization of the objective function values!
%
% Input:
% - F                   - A matrix of M x l, where M is the number
%                         of objectives, and l is the number of
%                         objective function value vectors of the
%                         solutions.
% - offset_fraction     - Optional: The shift of the reference point
%                         relative to the difference between the 
%                         minimum and maximum function value for
%                         each direction (default: 0)
%
% Output:
% - ref_point           - The reference point, a M x 1 vector.
%
% Author: Johannes W. Kruisselbrink
% Last modified: March 17, 2011

	if (nargin < 2), offset_fraction = 0; end

	ref_point = max(F,[],2) + offset_fraction * (max(F,[],2) - min(F,[],2));

end
