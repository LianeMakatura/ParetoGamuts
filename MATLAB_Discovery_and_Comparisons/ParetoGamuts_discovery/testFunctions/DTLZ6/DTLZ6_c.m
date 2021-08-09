function [fun2, deriv2, f2] = DTLZ6_c(N_Order)

A = sym('x', [1 N_Order]);
fun2 = @DTLZ6_c_helper;
f2(A) = DTLZ6_c_helper(A.');
A_cell = num2cell(A);
deriv2(A) = gradient(f2(A_cell{:}), A);
