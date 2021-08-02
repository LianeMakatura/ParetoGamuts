function [ out ] = DTLZ2a_helper( x )
    [X,Y] = size(x);

    x_first = x(1,:);
    x_second = x(2,:);
     
    summation = zeros(1,Y);
    
    for i=3:X
        xi = x(i,:);
        term = (xi-0.5).^2 - cos(20*pi*(xi-0.5));
        summation = summation + term;
    end

    g = 100*(X - 2 + summation);
   
    out = (0.5 * x_first.*x_second.*(1+g))/(120*(X-2));
end

