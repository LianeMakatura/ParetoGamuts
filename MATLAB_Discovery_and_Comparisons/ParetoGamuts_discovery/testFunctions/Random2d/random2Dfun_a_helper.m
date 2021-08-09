function out = random2Dfun_a_helper(x)
scaledX = x*2*pi;
[X,Y] = size(x);

coefficients = rand(1, X);
offsets = rand(1,X);


out = zeros(1,Y);

for i = 1:Y
    
    term = coefficients(i)* sin(scaledX(:,i) + offsets(i));
    
    out = out + term;
end

end

