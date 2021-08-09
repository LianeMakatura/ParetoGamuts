function [ out ] = DTLZ1b_helper( x )
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
   
    out = cos(x_first*(pi/2)).*sin(x_second*pi/2).*(1+g);
    out = out/(220*(X-2));
    
end

