% L0325002_Arya Lani Riyonggani Nurdian
clc; clear;
format long g;

% Matriks A dan Vektor b
A = [
     5, -1, -1,  0,  0;
    -1,  5, -1, -1,  0;
    -1, -1,  4, -1, -1;
     0,  0,  1,  4, -2;
     0,  1, -1,  1,  4
];

b = [-1; 2; 6; 2; -1];

max_iter = 20; % Petunjuk Soal 1
tol = 0.001;   % Petunjuk Soal 1

fprintf('=== PENYELESAIAN SPL (SOAL 1) ===\n');

% 1. Metode Gauss
[x_gauss, time_gauss] = EliminasiGauss(A, b);
fprintf('Hasil Eliminasi Gauss:\n'); disp(x_gauss');

% 2. Metode Gauss-Jordan
[x_gj, time_gj] = EliminasiGaussJordan(A, b);
fprintf('Hasil Eliminasi Gauss-Jordan:\n'); disp(x_gj');

% 3. Dekomposisi LU
[x_lu, time_lu] = LU_Solusi(A, b);

% Menampilkan nilai pengali
[L_lu, U_lu, P_lu] = LU_fwdelim(A);
fprintf('Nilai pengali (m):\n');
for i = 1 : size(A, 1) - 1
    for h = i+1 : size(A, 1)
        fprintf('m%d%d = %g\n', h, i, L_lu(h, i));
    end
end

fprintf('Hasil Dekomposisi LU:\n'); disp(x_lu');

% 4. Metode Jacobi
[x_jacobi, time_jacobi] = Jacobi(A, b, tol, max_iter);
fprintf('Hasil Iterasi Jacobi:\n'); disp(x_jacobi');

% 5. Metode Gauss-Seidel
[x_gs, time_gs] = IterasiGaussSeidel(A, b, tol, max_iter);
fprintf('Hasil Iterasi Gauss-Seidel:\n'); disp(x_gs');


fprintf('\n=== PERBANDINGAN GALAT/ERROR (SOAL 2) ===\n');
% Menggunakan Eliminasi Gauss sebagai nilai eksak pembanding
nilai_eksak = x_gauss;

galat_jacobi = abs(nilai_eksak - x_jacobi);
galat_gs = abs(nilai_eksak - x_gs);

fprintf('Galat Abstrak Jacobi (abs(Eksak - Hampiran):\n');
disp(galat_jacobi');

fprintf('Galat Abstrak Gauss-Seidel (abs(Eksak - Hampiran)):\n');
disp(galat_gs');


fprintf('\n=== PERBANDINGAN RUNNING TIME (SOAL 3) ===\n');
fprintf('---------------------------------------------------\n');
fprintf('Metode                | Running Time (detik)\n');
fprintf('---------------------------------------------------\n');
fprintf('Eliminasi Gauss       | %.6f\n', time_gauss);
fprintf('Eliminasi Gauss-Jordan| %.6f\n', time_gj);
fprintf('Dekomposisi LU        | %.6f\n', time_lu);
fprintf('Iterasi Jacobi        | %.6f\n', time_jacobi);
fprintf('Iterasi Gauss-Seidel  | %.6f\n', time_gs);
fprintf('---------------------------------------------------\n');
