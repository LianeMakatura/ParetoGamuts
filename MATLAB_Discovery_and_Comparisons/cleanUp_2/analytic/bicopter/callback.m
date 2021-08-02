function [F, J, H] = callback(x, N, fm, dfm)
    ep = 1e-6;
    [q, u, dt, m, r] = unpack(x, N);
    
    F = fm(q, u, dt, m, r);
    J = dfm(q, u, dt, m, r);
    
    H = zeros(41, 41);
    
    for i = 1:1:41
        x_high = x;
        x_high(i) = x_high(i) + ep;
        x_low = x;
        x_low(i) = x_low(i) - ep;
        
        [q_high, u_high, dt_high, m_high, r_high] = unpack(x_high, N);
        [q_low, u_low, dt_low, m_low, r_low] = unpack(x_low, N);

        
        df_high = dfm(q_high, u_high, dt_high, m_high, r_high);
        df_low = dfm(q_low, u_low, dt_low, m_low, r_low);
        
        H(:, i) = (df_high - df_low) / (2 * ep);               
    end        
end