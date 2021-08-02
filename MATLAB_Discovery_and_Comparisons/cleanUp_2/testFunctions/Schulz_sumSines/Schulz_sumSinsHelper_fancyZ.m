
function output = Schulz_sumSinsHelper_fancyZ(randA, randB,  origMin, origMax, x, z)

nSamples = size(x,2);
rD = size(x,1);

if exist('z', 'var') %% we should fix z (overwrite x(rD)
    x(rD, :) = z;
end

% out = zeros(nSamples, 1);
output = [];

% disp("FUNCTION CALL");
z_val = x(rD,:);

for j=1:nSamples
    res = 0;
    for i=1:rD % this should be RD

        iters =3*(1:size(randA,2));

        %sum sign for one design parameter
        y = x(i, j).'; % a row of the ith design variables  
        a2 = randA(i,:);
        b2 = randB(i,:);

        v1 = iters.'*y.';

        v2 = v1 + repmat(b2.^(2*z_val(j)).', 1, size(y,1));
        
        sumSin = (a2*cos(v2));
        res = res + sumSin;
    end
    
    output = [output; res];

%     for i=1:rD % this should be RD
% 
%         iters =3*(1:size(randA,2));
% 
%         %sum sign for one design parameter
%         y = x(i, :).'; % a row of the ith design variables  
%         a2 = randA(i,:);
%         b2 = randB(i,:);
% 
% 
%         v1 = iters.'*y.';
%         v2 = v1 + repmat(b2.^z_val(1).', 1, size(y,1));
%         %v2 = bsxfun(@plus,v1 , b2');
%         sumSin = (a2*cos(v2));
%         out = out + sumSin';
%     end
end

output = (output - origMin)/(origMax - origMin);

end
