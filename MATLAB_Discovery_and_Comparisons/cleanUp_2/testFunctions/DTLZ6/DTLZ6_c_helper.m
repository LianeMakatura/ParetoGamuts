function [ out ] = DTLZ6_c_helper( x )
    [X,Y] = size(x);
     
    summation = zeros(1,Y);
    for i=3:X
        summation = summation +  x(i,:);
    end
    
    g = 1 + (9/(X-2))*summation;
    
    summation_2 = zeros(1,Y);
    for i = 1:2
        term_1 = x(i,:)./(1 + g);
        term_2 = 1 + sin(3*pi*x(i,:));
        fullTerm = term_1.*term_2;
        summation_2 = summation_2 + fullTerm;
    end
   
    h = 3 - summation_2;
    
    out = (1 + g).*h;
    out = out /30;
    
end

