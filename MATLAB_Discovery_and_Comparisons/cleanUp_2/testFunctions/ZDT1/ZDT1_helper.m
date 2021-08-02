function [ out ] = ZDT1_helper( x )

    x_first = x(1,:);
    [X,Y] = size(x);
    summation = zeros(1,Y);
    
    for i=2:X
        xi = x(i,:);
        summation = summation + xi;
    end
    
    g = 1 + (9*summation/(X-1));
    
    h = 1-sqrt(x_first./g);
    out = g.*h/10;
end

