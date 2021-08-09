function [fun2, deriv2, f2] = random2Dfun_a( N_Order)

A = sym('x', [1 N_Order]);
fun2 = @random2Dfun_a_helper;
f2(A) = random2Dfun_a_helper(A.');
A_cell = num2cell(A);

deriv2(A) = gradient(f2(A_cell{:}), A);

end