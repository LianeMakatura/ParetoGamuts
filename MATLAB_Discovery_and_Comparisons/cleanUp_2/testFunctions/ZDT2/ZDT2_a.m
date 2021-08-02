function [fun2, deriv2, f2] = ZDT2_a( N_Order)

A = sym('x', [1 N_Order]);
fun2 = @ZDT2_helper;
f2(A) = ZDT2_helper(A.');
A_cell = num2cell(A);
deriv2(A) = gradient(f2(A_cell{:}), A);


% deriv2
%fun = @(x,y) (0.*x + 0.*y + 0);
