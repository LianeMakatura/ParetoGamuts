function [fun2, deriv2, f2] = Schulz_ZDT4_a(N_Order, z)

A = sym('x', [1 N_Order]);
fun2 = @Schulz_ZDT4_helper;
f2(A) = Schulz_ZDT4_helper(A.', z);
A_cell = num2cell(A);
deriv2(A) = gradient(f2(A_cell{:}), A);

% deriv2
%fun = @(x,y) (0.*x + 0.*y + 0);
