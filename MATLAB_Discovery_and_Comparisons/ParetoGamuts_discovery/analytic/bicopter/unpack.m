function [q, u, dt, m, r] = unpack(x, N)

    q = x(1:6);
    u = x(7:7 + 2*N - 1);
    dt = x(8 + 2*N - 1);
    m = x(9 + 2*N - 1);
    r = x(10 + 2*N - 1);
end