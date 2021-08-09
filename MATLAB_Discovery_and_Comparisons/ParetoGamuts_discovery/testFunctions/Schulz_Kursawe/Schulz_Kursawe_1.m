function [fun, deriv2, f] = Schulz_Kursawe_1( N_Order, z)

A = sym('x', [1 N_Order]);

if exist('z', 'var')        % constrained z
    fun = @(x,z)Schulz_Kursawe_1_helper(x,z);
    f(A) = Schulz_Kursawe_1_helper(A.', z);
else                        % unconstrained z
    fun = @(x)Schulz_Kursawe_1_helper(x);
    f(A) = Schulz_Kursawe_1_helper(A.');
end

A_cell = num2cell(A);
deriv2(A) = gradient(f(A_cell{:}), A);

end

% deriv2

%fun = @(x,y) (0.*x + 0.*y + 0);
