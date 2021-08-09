function [des_extents] = appGetExtents(centerPt, explorationDirections, ...
                                        mFunc, rD, rd, rA, userParams)    
    % figure out how far to push based on the jacobian
    xVals = centerPt(1, 1:rD);
    zVals = centerPt(1, rD+1:rD+rA);
    J = mFunc.evalJacobian('both', xVals, zVals, userParams.fd_eps, userParams.fd_type);

    mvmt = J*explorationDirections;
    
    % go for unit length movement: normalize movement then use the same
    % factor on the exploration directions
    mvmt_norm = sqrt(sum(mvmt.^2,1)); % norm of each column
    
    % in direction i, go from -des_extents(i):des_extents(i)
    des_extents = (1./mvmt_norm)*5;
end

