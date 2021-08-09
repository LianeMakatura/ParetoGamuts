function [ out ] = ZDT6_a_helper( x )

    x_first = x(1,:);
    first_term = exp(-4*x_first);
    second_term = power(sin(6*pi*x_first),6);
    out = 1 - first_term .* second_term;

end
