function [fun, deriv, f] = Schulz_appVarDim(N_Order)
    % SCHULZ_APPVARDIM Summary of this function goes here
    %   Detailed explanation goes here
    % assume z is 1 dimensional, and contained in the last element of the
    % design var vector
    
    A = sym('x', [1 N_Order]);
    A_cell = num2cell(A);

    fun = @(A) (A(N_Order,:));  % return the app var -- last element of the design vector
    f(A) = A(N_Order);
    deriv(A) = gradient(f(A_cell{:}), A);
end

