function [ out ] = DTLZ3c_helper( x )
    [X,Y] = size(x);

    x_first = x(1,:);
     
    summation = zeros(1,Y);
    
    for i=3:X
        xi = x(i,:);
        term = (xi-0.5).^2 - cos(20*pi*(xi-0.5));
        summation = summation + term;
    end
    
    g = 100*(X - 2 + summation);
   
    out = sin(x_first*(pi/2)).*(1+g);
    out = out/(220*(X-2));

end

