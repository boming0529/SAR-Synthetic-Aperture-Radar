% MATLAB 高斯白噪音範例

% 設定隨機噪音的參數
mu = 0;        % 均值
sigma = 1;     % 標準差
num_samples = 1000;  % 樣本數量

% 生成高斯白噪音
noise = mu + sigma * randn(num_samples, 1);

% 繪製噪音
figure;
plot(noise, 'Color', [0.5 0.5 0.5]); % 灰色線條
title('Gaussian White Noise');
xlabel('Sample index');
ylabel('Amplitude');


% SAR 圖像的 2-D 高斯白噪音生成

% 設定參數
Fr = 60e6;        % Range sampling rate, Hz
Naz = 1024;       % Number of range lines (rows)
Nrg = 320;        % Samples per range line (columns)
Fa = 200;         % Azimuth sampling rate or PRF, Hz

% 均值和標準差
mu = 0;           % 均值
sigma = 1;        % 標準差

% 生成 2-D 高斯白噪音
noise = mu + sigma * randn(Naz, Nrg);

% 繪製噪音
figure;
imagesc(noise);
colormap(gray);   % 使用灰階色彩
colorbar;
title('2-D Gaussian White Noise for SAR Image');
xlabel('Range Samples');
ylabel('Azimuth Lines');


% 三維繪圖
figure;
surf(noise, 'EdgeColor', 'none'); % 去除邊緣線以便更清晰地看到色彩
colormap(jet);  % 使用 jet 色圖增強視覺效果
colorbar;       % 顯示色條
caxis([-3 3]);  % 調整色彩軸範圍以使色條更明顯
xlabel('Range Samples');
ylabel('Azimuth Lines');
zlabel('Amplitude');
title('3D Visualization of Gaussian White Noise for SAR Image');