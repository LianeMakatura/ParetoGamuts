function [ out ] = Schulz_Kursawe_2_helper( x, z )
% x=x(1:3,:);


    [X,Y] = size(x);
    x(:,:) = 10*x(:,:) - 5; % rescale design vars to [0,1]
    rD = size(x,1);

    if exist('z', 'var') %% we should fix z (overwrite x(rD))
        x(rD, :) = z;
    end
    
    summation = zeros(1,Y);
    
    for i=1:X
        xi = x(i,:);
        additive = abs(xi).^0.8 + (5 * sin(xi.^3));
        
        summation = summation + additive;
    end
   
    out = summation;
    out = out./37.0;
    out = out + (12.0/37.0);
end

