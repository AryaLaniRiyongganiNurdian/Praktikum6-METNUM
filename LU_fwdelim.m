% L0325002_Arya Lani Riyonggani Nurdian
function [L, U, P] = LU_fwdelim(A)
[n, l] = size(A);
L = zeros(n);
P = eye(n);
for i = 1 : n-1
    [pivot, k] = max(abs(A(i:n, i)));
    k = k + i - 1;
    if abs(A(i, i)) < 1e-8 * pivot
        A([i k], :) = A([k i], :);
        P([i k], :) = P([k i], :);
        L([i k], 1:i-1) = L([k i], 1:i-1);
    end
    if abs(A(i, i)) < 1e-12
        error('Matriks singular, solusi tidak tunggal');
    end
    for h = i+1 : n
        m = A(h, i) / A(i, i);
        L(h, i) = m;
        A(h, :) = A(h, :) - m * A(i, :);
    end
end
if abs(A(n, n)) < 1e-12
    error('Matriks singular, solusi tidak tunggal');
end
L = L + eye(n);
U = A;
end
