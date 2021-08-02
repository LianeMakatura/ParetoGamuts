function [fun, deriv2, f] = Kursawe_2( N_Order)

A = sym('x', [1 N_Order]);
fun = @Kursawe_2_helper;
f(A) = Kursawe_2_helper(A.');
A_cell = num2cell(A);
deriv2(A) = gradient(f(A_cell{:}), A);
