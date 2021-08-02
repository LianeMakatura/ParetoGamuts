
function out = Schulz_sumSinsHelper(randA, randB,  origMin, origMax, x, z)

nSamples = size(x,2);
rD = size(x,1);

if ~exist('z', 'var') %% we're in the unconstrained z case, where z is a design var
    rD = rD-1;
    x = x(1:rD); % consider first two design vars
    z = x(rD);
    disp('sumSinsHelper unconstrained z!')
end

out = zeros(nSamples, 1);

for i=1:rD % this should be RD
    
    iters = 3*(1:size(randA,2));
    %sum sign for one design parameter
    y = x(i, :).'; % a row of the ith desing variables  
    a2 = randA(i,:);
    b2 = randB(i,:);

    v1 = iters.'*y.';

    v2 = v1 + repmat(b2.', 1, size(y,1));
    %v2 = bsxfun(@plus,v1 , b2');
    sumSin = (a2*cos(v2));
    out = out + sumSin';
   
    

end

 out = (out - origMin)/(origMax - origMin) + z; % +z was not there in original
end
