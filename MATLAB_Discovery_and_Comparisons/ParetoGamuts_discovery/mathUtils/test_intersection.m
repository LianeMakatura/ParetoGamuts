n1 = [1,0,0];
n2 = [0,0,1];

p0 = [0,0,0];
p1 = [0,0,0]; % pt on plane 1
p2 = [4,0,2]; % pt on plane 2

closest_point = plane_intersection(n1, p1, n2, p2, p0);

disp(closest_point)