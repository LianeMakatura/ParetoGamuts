function [closest_point] = plane_intersection(n1, p1, n2, p2, p0)
    % n1 is normal of first plane
    % p1 pt on first plane
    % n2 normal of second plane
    % p2 pt on second plane
    % p0 is target pt
    
    % closest point -- point along intersection of the 2 planes which is
    %                   closest to p0
    
    M = zeros(5);

    M(1,1) = 2;
    M(2,2) = 2;
    M(3,3) = 2;

    M(4, 1:3) = n1;
    M(5, 1:3) = n2;
    M(1:3, 4) = n1';
    M(1:3, 5) = n2';

    b = zeros(5, 1);
    b(1:3) = 2*p0';
    b(4) = dot(p1, n1);
    b(5) = dot(p2, n2);

    x = pinv(M)*b;

    closest_point = x(1:3);
end
