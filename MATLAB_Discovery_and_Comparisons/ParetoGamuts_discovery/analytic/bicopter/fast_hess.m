function [hess] = fast_hess(metric, jacob, design_variables, app_variables)
    eps = 1e-6;
    F = metric(design_variables, app_variables);
    J = jacob(design_variables, app_variables);
    x = [design_variables; app_variables]; % 33 x 1
    [rD, ~] = size(design_variables);
    [rA, ~] = size(app_variables);
    numVars = rD + rA;
    hess = [];
    for i = 1:numVars
        x_high = x;
        x_high(i) = x_high(i) + eps;
        x_low = x;
        x_low(i) = x_low(i) - eps;
        df_high = jacob(x_high(1:rD), x_high(rD+1:numVars)); % 1 x 33
        df_low = jacob(x_low(1:rD), x_low(rD+1:numVars));
        hess = [hess; (df_high - df_low) / (2 * eps)];
    end
end
