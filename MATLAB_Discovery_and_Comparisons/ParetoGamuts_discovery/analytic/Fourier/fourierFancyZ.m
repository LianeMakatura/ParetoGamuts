function [f]= fourierFancyZ(N_Order, rD, rA)
% z always a design var in this case; but if it's passed in separately, 
% x(rD)=z is fixed at that value 
dim = rD+rA;
origMin = -dim;
origMax = dim;

randA = rand(dim, N_Order);
randA = randA./ sum(randA,2); 
randB = pi*rand(dim, N_Order);

totalMin = 0;
totalMax = 0;

f_0 = 0;

for i=1:dim
    iters =3*(1:size(randA,2));
    a2 = randA(i,:);
    b2 = randB(i,:);   
    g_min = @(x)(a2*cos(iters.'*x + b2.'));
    g_max = @(x)(-a2*cos(iters.'*x + b2.'));
    
    f_0 = f_0 + g_min(0);

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

xSymb = sym('x', [1 rD]);
zSymb = sym('z', [1 rA]);

f(xSymb, zSymb) = fourierFancyZ_helper(randA, randB, totalMin - 0.1, totalMax + 0.1, xSymb.', zSymb.');
% 
% deriv(A) = gradient(f, [xSymb, zSymb]);
% hess(A) = hessian(f, [xSymb, zSymb]);
% fun = @(x, z)fourierFancyZ_helper(randA, randB, totalMin - 0.1, totalMax + 0.1, x, z);

end

