function [ out ] = Schulz_ZDT2_helper( x,z )

    x_first = x(1,:);
    [X,Y] = size(x);
    summation = zeros(1,Y);
    
    for i=2:X-1
        xi = x(i,:);
        summation = summation + xi;
    end
    summation = summation + z; 
    
    g = 1 + (9*summation/(X-1));
    h = 1-(x_first./g).^2;
   
    out = g.*h/10;
%     out = g.*h/z;
end

