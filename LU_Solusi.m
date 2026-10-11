% L0325002_Arya Lani Riyonggani Nurdian
function [x, elapsed_time] = LU_Solusi(A, b)
    t_start = tic();
    [L, U, P] = LU_fwdelim(A);
    [n, m] = size(L);
    b = P * b;
    y = zeros(n, 1);
    x = zeros(n, 1);
    y(1) = b(1) / L(1, 1);
    for i = 2 : n
        y(i) = (b(i) - L(i, 1:i-1) * y(1:i-1)) / L(i, i);
    end
    x(n) = y(n) / U(n, n);
    for i = n-1 : -1 : 1
        x(i) = (y(i) - U(i, i+1:n) * x(i+1:n)) / U(i, i);
    end
    elapsed_time = toc(t_start);
end
