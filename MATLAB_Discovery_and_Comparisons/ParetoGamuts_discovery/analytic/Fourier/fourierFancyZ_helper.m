
function output = fourierFancyZ_helper(randA, randB,  origMin, origMax, x, z)

nSamples = size(x,2);
rD = size(x, 1);
rA = size(z, 1);
dim = rD + rA;

% out = zeros(nSamples, 1);
output = [];

z_val = z(1:rA,:); %% IMPLEMENTATION CURRENTLY ONLY FOR RA=1

for j=1:nSamples
    res = 0;
    for i=1:dim % this should be RD

        iters =3*(1:size(randA,2));

        %sum sign for one design parameter
        if i <= rD
            y = x(i, j).'; % a row of the ith design variables  
        else
            y = z(i-rD, j).'; % a row of the ith application variables
        end
        a2 = randA(i,:);
        b2 = randB(i,:);

        v1 = iters.'*y.';

        v2 = v1 + repmat(b2.^(2*z_val(j)).', 1, size(y,1));
        
        sumSin = (a2*cos(v2));
        res = res + sumSin;
    end
    
    output = [output; res];

end

output = (output - origMin)/(origMax - origMin);

end
