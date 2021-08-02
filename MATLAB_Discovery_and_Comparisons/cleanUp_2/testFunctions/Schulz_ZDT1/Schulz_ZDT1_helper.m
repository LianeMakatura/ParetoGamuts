function [ out ] = Schulz_ZDT1_helper( x , z)

    x_first = x(1,:);
    [X,Y] = size(x);
    summation = zeros(1,Y);
    
    for i=2:X-1
        xi = x(i,:);
        summation = summation + xi;
    end
    summation = summation + z; 
    
    g = 1 + (9*summation/(X-1));
    
    h = 1-sqrt(x_first./g);
    out = g.*h/10; % original line
%     out = g.*h/z;
end

