function [ out ] = ZDT3_helper( x )

    x_first = x(1,:);
    [X,Y] = size(x);
    summation = zeros(1,Y);
    
    for i=2:X
        xi = x(i,:);
        summation = summation + xi;
    end
    
%     G
    g = 1 + (9*summation/(X-1));

%     H
    first_part = 1-sqrt(x_first./g);
    div = (x_first./g);
    trig = sin(10*pi*x_first);
    
    second_part = div.*trig;
    h = first_part - second_part;
   
    out = (1 + g.*h)/11;
end

