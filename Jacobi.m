% L0325002_Arya Lani Riyonggani Nurdian
function [x, elapsed_time] = Jacobi(A, b, tol, max_iter)
    t_start = tic();
    [n, m] = size(A);

    for i = 1:n
        xlama(i) = b(i) / A(i, i);
    end
    xlama = xlama';

    C = -A;
    for i = 1:n
        C(i, i) = 0.0;
        C(i, :) = C(i, :) / A(i, i);
        d(i, 1) = b(i) / A(i, i);
    end

    i = 1;
    while (i <= max_iter)
        xbaru = C * xlama + d;
        if (abs(xbaru - xlama) <= tol)
            x = xbaru;
            elapsed_time = toc(t_start);
            disp('Jacobi method konverge');
            return;
        else
            xlama = xbaru;
        end
        i = i + 1;
    end

    x = xbaru;
    elapsed_time = toc(t_start);
    disp('Jacobi method not konverge within max_iter');
end
