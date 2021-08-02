function [fun, deriv2, f] = Kursawe_1( N_Order)

A = sym('x', [1 N_Order]);
fun = @Kursawe_1_helper;
f(A) = Kursawe_1_helper(A.');
A_cell = num2cell(A);
deriv2(A) = gradient(f(A_cell{:}), A);

% deriv2

%fun = @(x,y) (0.*x + 0.*y + 0);

