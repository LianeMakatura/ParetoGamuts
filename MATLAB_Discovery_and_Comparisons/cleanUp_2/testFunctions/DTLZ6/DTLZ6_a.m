function [fun, deriv, f] = DTLZ6_a(N_Order)

A = sym('x', [1 N_Order]);
A_cell = num2cell(A);

fun = @(A) (A(1,:));
f(A) = A(1);
deriv(A) = gradient(f(A_cell{:}), A);