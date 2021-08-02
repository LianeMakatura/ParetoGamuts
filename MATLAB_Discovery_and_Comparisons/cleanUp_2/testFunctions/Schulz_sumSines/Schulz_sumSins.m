function [fun, deriv, f]= Schulz_sumSins(N_Order, rD, origMin, origMax, z)

if ~exist('z', 'var') %% we're in the unconstrained z case, where z is a design var
    rD = rD-1;
    disp('sumSins unconstrained z!')
end

% a = rand(1, No); 

randA = rand(rD, N_Order);
randA = randA./ sum(randA,2); 
randB = pi*rand(rD, N_Order);

disp('randA')
disp(randA)
disp('randB')
disp(randB)

totalMin = 0;
totalMax = 0;

f_0 = 0;
for i=1:rD % this should be RD
    
    iters =3*(1:size(randA,2));
    a2 = randA(i,:);
    b2 = randB(i,:);   
    g_min = @(x)(a2*cos(iters.'*x + b2.'));
    g_max = @(x)(-a2*cos(iters.'*x + b2.'));
    
    
    f_0 = f_0 + g_min(0);
    
 %   h_min = @(x)(-a2.*iters*sin(iters.'*x + b2.'));
 %   h_max = @(x)(a2.*iters*sin(iters.'*x + b2.'));
    
    
    minVal = inf;
    maxVal = -inf;
    options = optimoptions('fmincon','Display', 'off');
    for x0=0:0.01:1
        [~,lmin]  = fmincon(g_min,x0,[],[],[],[], 0,1, [], options);
        [~,lmax]  = fmincon(g_max,x0,[],[],[],[], 0,1, [], options);
        if (lmin < minVal)
            minVal = lmin;
        end
        if (-lmax > maxVal)
            maxVal = -lmax;
        end
    end
    totalMin = totalMin + minVal;
    totalMax = totalMax + maxVal;
end

if exist('z', 'var')
    A = sym('x', [1 rD]);

    fun = @(x, z)Schulz_sumSinsHelper(randA, randB, totalMin - 0.1, totalMax + 0.1, x, z);
    f(A) = Schulz_sumSinsHelper(randA, randB, totalMin - 0.1, totalMax + 0.1, A.', z);
else % unconstrained z
    A = sym('x', [1 rD+1]);

    fun = @(x)Schulz_sumSinsHelper(randA, randB, totalMin - 0.1, totalMax + 0.1, x);
    f(A) = Schulz_sumSinsHelper(randA, randB, totalMin - 0.1, totalMax + 0.1, A.');
end

deriv(A) = gradient(f, A);

end
% i =3:N_Order;
% No = size(i,2);
% a = rand(1, No); 
% D = 5;
% b = D*rand(1, No) - D/2;
% a2 = rand(1, No);
% b2 = D*rand(1, No) - D/2;
% a3 = rand(1, No);
% b3 = D*rand(1, No) - D/2;
% 
% 
% 
% 
% fun = @(x) (a*sin( bsxfun(@plus,i'*x(:,1) , b'))).*(a2*cos( bsxfun(@plus,i'*y , b2'))).*(a3*cos( bsxfun(@plus,i'*z , b3')));
% 
% 
% syms sx sy sz
% f(sx, sy, sz) = (a*sin(i'*sx + b')).*(a2*cos( i'*sy + b2')).*(a3*cos( i'*sz + b3'));
% deriv(sx, sy, sz)  = gradient(f(sx,sy, sz), [sx, sy, sz]);
% 
% 
% 


