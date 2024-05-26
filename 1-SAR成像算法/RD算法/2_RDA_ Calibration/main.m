clear;
clc;

% Add the path to the directory containing calculate_image_power function
relativePath = '../2_小斜視角_RDA_White_noise';
addpath(relativePath);

load('corss_section.mat', 'corrected_img');

[max_col_value, max_col_index] = max(abs(corrected_img));           
[max_rol_value, max_rol_index] = max(max(abs(corrected_img)));
row_max = max_col_index(max_rol_index);
column_max = max_rol_index;

original_size = 15;
scale_factor = 7;
NN = floor((scale_factor * original_size - 1) / 2);
DN = corrected_img(row_max - NN: row_max + NN, column_max - NN: column_max + NN); 

Fr = 60e6;                  % 距離採樣率   (Range sampling rate, MHz)
Fa = 200;                   % 方位採樣率   (Azimuth sampling rate or PRF, HZ)
Vr = 150;                   % 雷達有效速度 (Effective rader velocity, m/s)
c = 3e8;                    % 光速(Speed of light)
delta_a = 1 / Fa * Vr;
delta_b = 1 / Fr * c / 2;

epsilon_p = calculate_image_power(abs(DN), delta_a, delta_b);
disp(epsilon_p);


rmpath(relativePath);