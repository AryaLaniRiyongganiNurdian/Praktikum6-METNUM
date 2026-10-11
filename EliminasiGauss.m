% L0325002_Arya Lani Riyonggani Nurdian
function [x, elapsed_time] = EliminasiGauss(A, b)
    t_start = tic();
    [n, l] = size(A);
    for i = 1 : n-1
        [pivot, k] = max(abs(A(i:n, i)));
        k = k + i - 1;
        if abs(A(i, i)) < 1e-8 * pivot
            temp1 = A(i, :);
            temp2 = b(i, :);
            A(i, :) = A(k, :);
            b(i, :) = b(k, :);
            A(k, :) = temp1;
            b(k, :) = temp2;
        end
        if abs(A(i, i)) < 1e-12
            error('Matriks singular, solusi tidak tunggal');
        end
        for h = i+1 : n
            m = A(h, i) / A(i, i);
            A(h, :) = A(h, :) - m * A(i, :);
            b(h, :) = b(h, :) - m * b(i, :);
        end
    end
    if abs(A(n, n)) < 1e-12
        error('Matriks singular, solusi tidak tunggal');
    end
    x(n, :) = b(n, :) / A(n, n);
    for i = n-1 : -1 : 1
        x(i, :) = (b(i, :) - A(i, i+1:n) * x(i+1:n, :)) / A(i, i);
    end
    elapsed_time = toc(t_start);
end
