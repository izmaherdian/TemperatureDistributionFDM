% =========================================================================
% IK5025 - Analisis Variasi Material Logam (Fluks Kalor)
% Membuktikan bahwa k tidak mengubah profil T, tetapi mengubah Magnitudo q''
% =========================================================================

clear; clc; close all;

%% 1. PROSES SOLVER (Diringkas dari script sebelumnya)
L = 1.0; dx = 0.001; dy = dx;
Nx = round(L/dx) + 1; Ny = round(L/dy) + 1;
i1 = round(0.375/dx) + 1; i2 = round(0.625/dx) + 1;
j1 = round(0.375/dy) + 1; j2 = round(0.625/dy) + 1;

T = ones(Nx, Ny) * 47.5; 
T(1, :) = 100; T(Nx, :) = 20; T(:, 1) = 50; T(:, Ny) = 20;

tol = 1e-5; err = 1;
while err > tol
    T_old = T;
    for i = 2:Nx-1
        for j = 2:Ny-1
            if i > i1 && i < i2 && j > j1 && j < j2
                continue; 
            elseif i == i1 && j == j2
                T(i,j) = (2*T(i-1,j) + 2*T(i,j+1) + T(i+1,j) + T(i,j-1)) / 6;
            elseif i == i1 && j == j1
                T(i,j) = (2*T(i-1,j) + 2*T(i,j-1) + T(i+1,j) + T(i,j+1)) / 6;
            elseif i == i2 && j == j2
                T(i,j) = (2*T(i+1,j) + 2*T(i,j+1) + T(i-1,j) + T(i,j-1)) / 6;
            elseif i == i2 && j == j1
                T(i,j) = (2*T(i+1,j) + 2*T(i,j-1) + T(i-1,j) + T(i,j+1)) / 6;
            elseif i == i1 && j > j1 && j < j2
                T(i,j) = (2*T(i-1,j) + T(i,j+1) + T(i,j-1)) / 4;
            elseif i == i2 && j > j1 && j < j2
                T(i,j) = (2*T(i+1,j) + T(i,j+1) + T(i,j-1)) / 4;
            elseif j == j1 && i > i1 && i < i2
                T(i,j) = (2*T(i,j-1) + T(i+1,j) + T(i-1,j)) / 4;
            elseif j == j2 && i > i1 && i < i2
                T(i,j) = (2*T(i,j+1) + T(i+1,j) + T(i-1,j)) / 4;
            else
                T(i,j) = (T(i+1,j) + T(i-1,j) + T(i,j+1) + T(i,j-1)) / 4;
            end
        end
    end
    err = max(max(abs(T - T_old)));
end

% Hilangkan data di dalam isolator untuk proses gradien
T(i1+1:i2-1, j1+1:j2-1) = NaN;

%% 2. MENGHITUNG GRADIEN TEMPERATUR DAN FLUKS KALOR (HUKUM FOURIER)
% Menghitung gradien dT/dx dan dT/dy menggunakan fungsi built-in MATLAB
[dTdy, dTdx] = gradient(T, dy, dx);

% Properti Material (W/m.K)
k_Cu = 401; % Tembaga
k_CS = 43;  % Baja Karbon

% Fluks Kalor Vektor (q'' = -k * gradien T)
qx_Cu = -k_Cu * dTdx; qy_Cu = -k_Cu * dTdy;
qx_CS = -k_CS * dTdx; qy_CS = -k_CS * dTdy;

% Magnitudo Fluks Kalor Absolut |q''| = sqrt(qx^2 + qy^2)
q_mag_Cu = sqrt(qx_Cu.^2 + qy_Cu.^2);
q_mag_CS = sqrt(qx_CS.^2 + qy_CS.^2);

%% 3. PLOTTING PERBANDINGAN FLUKS KALOR
x_lin = linspace(0, L, Nx); y_lin = linspace(0, L, Ny);
[X, Y] = meshgrid(x_lin, y_lin);

% Membuat figure dengan ukuran spesifik dalam centimeter
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 28, 12]);

% --- SUBPLOT 1: BAJA KARBON ---
subplot(1,2,1);
contourf(X, Y, q_mag_CS', 25, 'LineColor', 'none');
colormap('hot'); 
c1 = colorbar; c1.Label.String = 'Fluks Kalor, q'''' (W/m^2)';
axis equal; axis tight;
title('Magnitudo Fluks Kalor: Baja Karbon (k = 43)', 'FontSize', 12);
xlabel('Jarak X (m)'); ylabel('Jarak Y (m)');
hold on;
rectangle('Position', [0.375, 0.375, 0.25, 0.25], 'FaceColor', 'k');
text(0.5, 0.5, 'Insulator', 'Color', 'w', 'HorizontalAlignment', 'center');
hold off;

% --- SUBPLOT 2: TEMBAGA ---
subplot(1,2,2);
contourf(X, Y, q_mag_Cu', 25, 'LineColor', 'none');
colormap('hot'); 
c2 = colorbar; c2.Label.String = 'Fluks Kalor, q'''' (W/m^2)';
axis equal; axis tight;
title('Magnitudo Fluks Kalor: Tembaga (k = 401)', 'FontSize', 12);
xlabel('Jarak X (m)'); ylabel('Jarak Y (m)');
hold on;
rectangle('Position', [0.375, 0.375, 0.25, 0.25], 'FaceColor', 'k');
text(0.5, 0.5, 'Insulator', 'Color', 'w', 'HorizontalAlignment', 'center');
hold off;

% --- PENGATURAN KERTAS UNTUK EXPORT AGAR TIDAK TERPOTONG ---
set(fig, 'PaperUnits', 'centimeters');
set(fig, 'PaperSize', [28 12]);         % Memaksa kertas menjadi landscape 28x12 cm
set(fig, 'PaperPosition', [0 0 28 12]); % Memenuhi seluruh area kertas

fprintf('Analisis Selesai');