% test finite differencing code on the lamp
mFunc = AppMappingFunction(29);
pDes = rand(1,mFunc.rD);
pApp = rand(1, mFunc.rA);
eps = 0.001;
fd_type = 'central';

disp('------- Testing fd jacobian')
wrt_vars = 'both';

fd_jac = fd_jacobian(mFunc, wrt_vars, pDes, pApp, eps, fd_type);

desCell = num2cell(pDes);
appCell = num2cell(pApp);
J = eval(mFunc.deriv(desCell{:}, appCell{:}));

figure;
imagesc(abs(fd_jac - J'));
colorbar
title('Absolute error: FD jacobian vs analytic')



disp('------- Testing fd hessian')
wrt_vars1 = 'both';
wrt_vars2 = 'both';

fd_hess = fd_hessian(mFunc, wrt_vars1, wrt_vars2, pDes, pApp, eps, fd_type);


desCell = num2cell(pDes);
appCell = num2cell(pApp);

h = mFunc.hess;
H = zeros(mFunc.rD+mFunc.rA,mFunc.rD+mFunc.rA,mFunc.rd);

for i = 1:mFunc.rd
    h_i = h{i};
    h_i_eval = double(h_i(desCell{:}, appCell{:}));
    H(:,:,i) = h_i_eval;
    
    figure;
    imagesc(abs(fd_hess{i} - h_i_eval));
    colorbar
    title('Absolute error: FD hessian vs analytic')
end

