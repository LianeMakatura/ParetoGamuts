function [projPatch] = proj_patch2subspace(patchPts, subspace)
%PROJ_PATCH2SUBSPACE Project a set of points to a subspace
%   @param patchPts -   M x N array of points, where each pts is on a row (M points, each N dimensional)
%   @param subspace     n x N array of n directions, each N dimensional
    projPatch = patchPts * subspace';
end

