t = Turbine(pwd, [1, 2]);
p = t.performanceMetric2();
p = matlabFunction(p);

pts = 0:0.01:1;
numpts = length(pts);
values = zeros(numpts, numpts);

for i=1:numpts % corresponds to y axis
    ival = pts(i);
    for j = 1:numpts % corresponds to x axis
        jval = pts(j);
        values(i,j) = p(ival,jval);
    end
end

[X, Y] = meshgrid(pts, pts);

figure
surf(X,Y,values);
ylabel('radius');
xlabel('windspeed');
zlabel('power');

[maxval,maxidx] = max(values(:));
[yidx, xidx] = ind2sub([numpts, numpts], maxidx);
fprintf("Max value is %0.6f at index (%d, %d)\n", maxval, yidx, xidx);

[minval,minidx] = min(values(:));
[yidx, xidx] = ind2sub([numpts, numpts], minidx);
fprintf("Min value is %0.6f at index (%d, %d)\n", minval, yidx, xidx);