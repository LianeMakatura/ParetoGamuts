function [slope] = slope2D(p1, p2)
    slope = (p2(2) - p1(2)) / (p2(1) - p1(1));
end