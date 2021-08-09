function [ out ] = ZDT6_b_helper( x )
    x_first = x(1,:);
    first_term = exp(-4*x_first);
    second_term = power(sin(6*pi*x_first),6);
    f1 = 1 - first_term .* second_term;
    
    [X,Y] = size(x);
    summation = zeros(1,Y);
    
    for i=2:X
        xi = x(i,:);
        additive = xi;
        summation = summation + additive;
    end
    
    term = (summation/(X-1)).^(1/4);
    
    g = 1 + 9*term;

    out = (1-(f1./g).^2)/10;
end