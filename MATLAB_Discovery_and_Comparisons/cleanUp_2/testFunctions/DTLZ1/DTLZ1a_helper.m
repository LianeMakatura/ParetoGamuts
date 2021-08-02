function [ out ] = DTLZ1a_helper( x )
    [X,Y] = size(x);

    x_first = x(1,:);
    x_second = x(2,:);
     
    summation = zeros(1,Y);
    
    for i=3:X
        xi = x(i,:);
        term = (xi-0.5).^2;
        summation = summation + term;
    end
    
    g = (summation);
   
    out = cos(x_first*(pi/2)).*cos(x_second*pi/2).*(1+g);
    out = out/2.25;
end

