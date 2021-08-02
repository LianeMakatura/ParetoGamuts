function [ out ] = Schulz_ZDT4_helper( x , z)

    x_first = x(1,:);
    [X,Y] = size(x);
    
    %rescaling design variables
    x(2:X,:) = 10*x(2:X,:) -5; 
    z = 10.*z - 5;
    summation = zeros(1,Y);
    
    for i=2:X-1
        xi = x(i,:);
        additive = xi.^2 - 10*cos(4*pi*xi);
        summation = summation + additive;
    end
    summation = summation + z.^2 - 10*cos(4*pi*z);
    
    g = 1 + (10*(X-1)) + summation;
    h = 1-sqrt(x_first./g);
    
    scalingFactor = X;%*30; % original
    %scalingFactor = X*z;
    
    out = g.*h/scalingFactor;
end

