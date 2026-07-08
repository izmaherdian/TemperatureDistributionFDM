% =========================================================================
% IK5025 - Analisis Variasi Ukuran Kisi (Grid Independence Test)
% Mengukur Performa Komputasi & Stabilitas Nilai Fisis (Probe Point)
% =========================================================================

clear; clc; close all;

%% 1. INISIALISASI PARAMETER PENGUJIAN
L = 1.0;                             % Dimensi pelat (m)
dx_array = [0.1, 0.05, 0.01]; % Variasi grid: Kasar, Menengah, Halus, Sangat Halus
num_tests = length(dx_array);

% Pre-alokasi untuk menyimpan hasil performa
waktu_komputasi = zeros(1, num_tests);
total_iterasi   = zeros(1, num_tests);
total_node      = zeros(1, num_tests);
T_results       = cell(1, num_tests);

% --- TITIK TINJAU (PROBE POINT) ---
x_probe = 0.75; 
y_probe = 0.75;
T_probe = zeros(1, num_tests); 

tol = 1e-4; % Gunakan 1e-4 
max_iter = 100000; % Batas aman iterasi dinaikkan karena grid 0.001 butuh banyak iterasi

fprintf('=======================================================================\n');
fprintf('          ANALISIS KINERJA GRID (GRID INDEPENDENCE TEST)\n');
fprintf('=======================================================================\n');
fprintf('PERINGATAN: Grid 0.001 m memiliki 1.002.001 node. Proses ini mungkin\n');
fprintf('memakan waktu belasan hingga puluhan menit. Silakan tunggu...\n');
fprintf('-----------------------------------------------------------------------\n');

%% 2. LOOPING KOMPUTASI UNTUK SETIAP VARIASI GRID
for k = 1:num_tests
    dx = dx_array(k);
    dy = dx;
    Nx = round(L/dx) + 1; 
    Ny = round(L/dy) + 1;
    total_node(k) = Nx * Ny;
    
    % Indeks batas isolator
    i1 = round(0.375/dx) + 1; i2 = round(0.625/dx) + 1;
    j1 = round(0.375/dy) + 1; j2 = round(0.625/dy) + 1;

    % Inisialisasi
    T = ones(Nx, Ny) * 47.5; 
    T(1, :) = 100; T(Nx, :) = 20; T(:, 1) = 50; T(:, Ny) = 20;

    err = 1; iter = 0;
    
    tic; % Mulai hitung waktu
    while err > tol && iter < max_iter
        T_old = T;
        iter = iter + 1;
        
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
    waktu_komputasi(k) = toc; % Selesai hitung waktu
    total_iterasi(k) = iter;
    
    % --- EKSTRAKSI DATA PROBE POINT ---
    i_probe = round(x_probe/dx) + 1;
    j_probe = round(y_probe/dy) + 1;
    T_probe(k) = T(i_probe, j_probe);
    
    T(i1+1:i2-1, j1+1:j2-1) = NaN; 
    T_results{k} = T;
end

% Menghitung Relative Error Suhu Probe terhadap grid terhalus (0.001)
rel_error_T = abs(T_probe - T_probe(end)) ./ T_probe(end) * 100;

% Mencetak Tabel Hasil
fprintf('Grid (m) | Total Node | Iterasi | Waktu (s) | T di (0.75, 0.75) | Galat (%%)\n');
fprintf('-----------------------------------------------------------------------\n');
for k = 1:num_tests
    fprintf('%8.3f | %10d | %7d | %9.4f | %15.4f C | %7.3f %%\n', ...
        dx_array(k), total_node(k), total_iterasi(k), waktu_komputasi(k), T_probe(k), rel_error_T(k));
end
fprintf('=======================================================================\n');

%% 3. PLOTTING VISUAL (KONTUR)
% Diperlebar menjadi 40 cm agar 4 subplot tidak tergencet
fig1 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 40, 10]);
for k = 1:num_tests
    dx = dx_array(k);
    Nx = round(L/dx) + 1; Ny = round(L/dx) + 1;
    x_lin = linspace(0, L, Nx); y_lin = linspace(0, L, Ny);
    [X, Y] = meshgrid(x_lin, y_lin);
    
    subplot(1, 4, k);
    contourf(X, Y, T_results{k}', 30, 'LineColor', 'none'); colormap('jet');
    axis equal; axis tight;
    title(sprintf('\\Delta x = %.3f m', dx), 'FontSize', 11);
    xlabel('X (m)'); ylabel('Y (m)');
    
    hold on;
    rectangle('Position', [0.375, 0.375, 0.25, 0.25], 'FaceColor', 'k');
    plot(x_probe, y_probe, 'w+', 'MarkerSize', 10, 'LineWidth', 2); % Tandai Probe Point
    if k == 4
        text(0.5, 0.5, 'Insulator', 'Color', 'w', 'HorizontalAlignment', 'center', 'FontSize', 9);
    end
    hold off;
end
cb = colorbar; cb.Position = [0.93 0.15 0.01 0.75]; cb.Label.String = 'Temperatur (^\circC)';

% Export Contour
set(fig1, 'PaperUnits', 'centimeters', 'PaperSize', [40 10], 'PaperPosition', [0 0 40 10]);
exportgraphics(fig1, 'perbandingan_grid_kontur.png', 'Resolution', 300);

%% 4. PLOTTING KURVA GRID INDEPENDENCE
fig2 = figure('Color', 'w', 'Position', [100, 100, 500, 400]);
% Sumbu X diubah menjadi skala Logaritmik agar titik 100, 400, 10.000, dan 1.000.000 bisa terlihat semua
semilogx(total_node, T_probe, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
grid on;
title('Kurva Independensi Kisi (\textit{Grid Independence})', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Jumlah Total \textit{Node} (Skala Log)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Temperatur di Probe $(0.75, 0.75)$ ($^\circ$C)', 'Interpreter', 'latex', 'FontSize', 12);

fprintf('Analisis Selesai!');