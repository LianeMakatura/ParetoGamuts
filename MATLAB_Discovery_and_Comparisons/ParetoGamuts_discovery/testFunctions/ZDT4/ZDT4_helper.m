function [ out ] = ZDT4_helper( x )

    x_first = x(1,:);
    [X,Y] = size(x);
    
    %rescaling design variables
    x(2:X,:) = 10*x(2:X,:) -5; 
    summation = zeros(1,Y);
    
    for i=2:X
        xi = x(i,:);
        additive = xi.^2 - 10*cos(4*pi*xi);
        summation = summation + additive;
    end
    
    g = 1 + (10*(X-1)) + summation;
    h = 1-sqrt(x_first./g);
    
    scalingFactor = X*30;
    
    out = g.*h/scalingFactor;
end

