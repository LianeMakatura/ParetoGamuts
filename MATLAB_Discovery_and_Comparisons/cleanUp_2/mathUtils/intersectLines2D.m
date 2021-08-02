function [x, y] = intersectLines2D(p1, m1, p2, m2)
% Find intersection of two lines, given the slope and a point (not
% necessarily the intercept). Derived from the two-point form.
    p1x = p1(1); p1y = p1(2);
    p2x = p2(1); p2y = p2(2);
    
    x = (p2y - p1y + m1*p1x - m2*p2x) / (m1 - m2);
    y = m2*(x - p2x) + p2y;
end