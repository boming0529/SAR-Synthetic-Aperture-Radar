%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            小斜視角（3.5°）
%                CSA
%      CSA in the Low Squint Case
%              點目標模擬
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update Nov 26 2023. 

%%
% --------------------------------------------------------------------------------------------------------
% 清除workspace和figure以及Command Window
% --------------------------------------------------------------------------------------------------------
clear;
close all;
clc;

%%
% --------------------------------------------------------------------------------------------------------
% 定義參數
% --------------------------------------------------------------------------------------------------------
% Book Table 6.1 : C-Band AirBorn SAR Parameters Used in the Simulations
% --------------------------------------------------------------------------------------------------------
R_nc = 20e3;                % 景中心斜距   (Slant range of scene center, Km)
Vr = 150;                   % 雷達有效速度 (Effective rader velocity, m/s)
Tr = 2.5e-6;                % 發射脈衝時寬 (Transmitted pulse duration, \mu s)
Kr = 20e12;                 % 距離向調频率 (Range FM rate, Hz/sec)
f0 = 5.3e9;                 % 雷達工作頻率 (Radar center frequency, GHz)
BW_dop = 80;                % 多普勒带寬   (Doppler bandwidth, Hz)
Fr = 60e6;                  % 距離採樣率   (Range sampling rate, MHz)
Fa = 200;                   % 方位採樣率   (Azimuth sampling rate or PRF, HZ)
Naz = 1024;          	    % 方位採樣數,即數據矩陣列數(Number of range lines(column))——這裡修改為1024。
Nrg = 320;                  % 距離採樣數,即數據矩陣行數(Sample per range line(row)）
sita_r_c = (3.5*pi)/180;	% 波束斜视角,3.5度,這裡轉為弧度(Beam squint angel, rad)

% --------------------------------------------------------------------------------------------------------
% 定義參數
% --------------------------------------------------------------------------------------------------------
c = 3e8;                   % 光速(Speed of light)
R0 = R_nc*cos(sita_r_c);   % 與R_nc相對應的最近斜距,根據(4.17)計算.(Range of closest approach)
Nr = Tr*Fr;                % 線性調頻信號採樣點數,課本找不到 (◕‿◕)?
BW_range = Kr*Tr;          % 距離向帶寬.(表4.3).(Range bandwidth)
lamda = c/f0;              % 波長(Wavelength)
fnc = 2*Vr*sin(sita_r_c)/lamda;             %多普勒中心頻率,根據(4.33)計算.(Doppler centroid frequency)
La_real = 0.886*2*Vr*cos(sita_r_c)/BW_dop;  %方位向天線長度,根據(4.36)計算.(Antenna length)  
beta_bw = 0.886*lamda/La_real;              %雷達3dB波束
La = beta_bw*R0;        % 合成孔徑長度
a_sr = Fr / BW_range;   % 距離向過採樣率,根據(4.22)計算.(Range oversampling ratio) 
a_sa = Fa / BW_dop;     % 方位向過採樣率,尚未找到.(Azimuth oversampling factor)
Mamb = round(fnc/Fa);   % 多普勒模糊,根據(12.20)計算,與課本公式不同.(Doppler ambiguity)

NFFT_r = Nrg;           % 距離向FFT長度,即數據矩陣列數
NFFT_a = Naz;           % 方位向FFT長度,即數據矩陣列數

R_ref = R0;             % (SRA)参考目标选在场景中心，其最近斜距为 R_ref
fn_ref = fnc;        	% (SRA)参考目标的多普勒中心频率

% --------------------------------------------------------------------------------------------------------
% 定義座標軸
% --------------------------------------------------------------------------------------------------------
% Book Figure 6.2 : Positions of three targets used in the simulation
% --------------------------------------------------------------------------------------------------------
% 模擬目標位置
% 以距離向作為x軸正方向
% 以方位向作為y軸正方向
% --------------------------------------------------------------------------------------------------------
delta_R0 = 0;       % 將目標A的波束中心穿越时刻，定義為方位向時間零點
delta_R1 = 120; 	% 目標A和目標B的方位向距離差，120m
delta_R2 = 80;      % 目標B和目標C的距離向距离差，80m

% 目標A
x1 = R0;                            % 目標A的距離向距離
y1 = delta_R0 + x1*tan(sita_r_c);	% 目標A的的方位向距離 

% 目標B
x2 = x1;            % 目標B和目標A的距離向距離相同
y2 = y1 + delta_R1; % 目標B的方位向距離

% 目標C
x3 = x2 + delta_R2;                 % 目標C和目標B有距離向的距離差,為delta_R2
y3 = y2 + delta_R2*tan(sita_r_c);  	% 目標C的方位向距離

% 定義目標A,B,C位置向量矩陣
x_range = [x1,x2,x3];
y_azimuth = [y1,y2,y3];

% 計算三個目標各自的波束中心穿越时刻
nc_1 = (y1-x1*tan(sita_r_c))/Vr;    % 目標A的波束中心穿越时刻
nc_2 = (y2-x2*tan(sita_r_c))/Vr;    % 目標B的波束中心穿越时刻
nc_3 = (y3-x3*tan(sita_r_c))/Vr;    % 目標C的波束中心穿越时刻
nc_target = [nc_1,nc_2,nc_3];       % 定義目標A,B,C穿越时刻矩陣

%% 
% --------------------------------------------------------------------
% 定義tau和eta矩陣
% --------------------------------------------------------------------
% 距離列矩陣
tau = 2*R0/c + ( -Nrg/2 : (Nrg/2-1) )/Fr;  % 距離時間軸,課本找不到(◕‿◕)
% 生成距離矩陣
tau_mtx = ones(Naz,1)*tau;                 % 距離時間軸矩陣,矩陣大小：Naz*Nrg

% 方位列矩陣
eta = ( -Naz/2: Naz/2-1 )/Fa;              % 方位时间轴
% 生成方位矩陣
eta_mtx = eta.'*ones(1,Nrg);               % 方位時間轴矩陣,矩陣大小：Naz*Nrg

%% 
% --------------------------------------------------------------------
% 定義f_tau和f_eta矩陣
% --------------------------------------------------------------------
fr = ( -NFFT_r/2 : NFFT_r/2-1 )*( Fr/NFFT_r );                 % 距离频率轴
fr_mtx = ones(Naz,1)*fr;                                       % 距离频率轴矩阵，大小：Naz*Nrg

fa = fnc + fftshift( -NFFT_a/2 : NFFT_a/2-1 )*( Fa/NFFT_a );   % 方位频率轴
fa_mtx = fa.'*ones(1,Nrg);                                     % 方位频率轴矩阵，大小：Naz*Nrg

%% 
% --------------------------------------------------------------------------------------------------------
% 6.3.1 Raw data
% --------------------------------------------------------------------------------------------------------
s_echo = zeros(Naz,Nrg);    % s_echo用来存放生成的回波(Raw Data)結果

A0 = 1;                     % 目標回波幅度,都設置為1
for k = 1:3                 % 生成k個目標的原始回波數據
    R_n = sqrt( (x_range(k).*ones(Naz,Nrg)).^2 + (Vr.*eta_mtx-y_azimuth(k).*ones(Naz,Nrg)).^2 );% 目標k的瞬時斜距(6.2)
    w_range = ((abs(tau_mtx-2.*R_n./c)) <= ((Tr/2).*ones(Naz,Nrg)));                            % 距離向包络,即距離窗(Range evelope)
    % ===================================================================================================   
    % 利用合成孔徑長度,直接構造矩形窗(其實這里只是限制數據範圍,沒有真正加窗）(◕‿◕)?
    
    w_azimuth = (abs(eta - nc_target(k)) <= (La/2)/Vr);    % 方位向包络(Azimuth evelope) 
    w_azimuth = w_azimuth.'*ones(1,Nrg);                   % 生成Naz*Nrg的矩陣
    % =====================================================================     
    % 下式就是生成某一個點目標(目標k)的回波信號

    s_k = A0.*w_range.*w_azimuth.*exp(-(1j*4*pi*f0).*R_n./c).*exp((1j*pi*Kr).*(tau_mtx-2.*R_n./c).^2);
    % 經過幾次循環,生成幾個點目標的回波信號,相加即可

    s_echo = s_echo + s_k;
end

% --------------------------------------------------------------------------------------------------------
% Make Figure 1
% --------------------------------------------------------------------------------------------------------
% Book Figure 6.3 : Simulated radar signal (raw) data with three target
% --------------------------------------------------------------------------------------------------------
figure('Name','Simulated radar signal (raw) data with three target', 'NumberTitle','on'); 
subplot(2,2,1);
imagesc(real(s_echo));
title('(a) 實部 Real Part');
xlabel('Range time domain(Samples)');
ylabel('Azimuth time domain(Samples)');

subplot(2,2,2);
imagesc(imag(s_echo));
title('(b) 虚部 Imaginary Part');
xlabel('Range time domain(Samples)');
ylabel('Azimuth time domain(Samples)');

subplot(2,2,3);
imagesc(abs(s_echo));
title('(c) 強度 Magnitude');
xlabel('Range time domain(Samples)');
ylabel('Azimuth time domain(Samples)');

subplot(2,2,4);
imagesc(angle(s_echo));
title('(d) 相位 Phase');
xlabel('Range Time domain(Samples)');
ylabel('Azimuth time domain(Samples)');

figure;
subplot(2,2,1);
imagesc(abs(fft(s_echo,[],1)));
title('RD 频谱幅度');
subplot(2,2,2);
imagesc(angle(fft(s_echo,[],1)));
title('RD 频谱相位');
subplot(2,2,3);
imagesc(abs(fft2(s_echo)));
title('二维频谱幅度');
subplot(2,2,4);
imagesc(angle(fft2(s_echo)));
title('二维频谱相位');
% colormap(gray);

%%
% --------------------------------------------------------------------
% 变换到距离多普勒域，进行“补余RCMC”
% --------------------------------------------------------------------
s_rd = s_echo.*exp(-1j*2*pi*fnc.*(eta.'*ones(1,Nrg))); 	% 数据搬移

S_RD = fft(s_rd,NFFT_a,1);  % 进行方位向傅里叶变换，得到距离多普勒域频谱

D_fn_Vr = sqrt(1-lamda^2.*(fa.').^2./(4*Vr^2));     % 大斜视角下的徙动因子，列向量
D_fn_Vr_mtx = D_fn_Vr*ones(1,Nrg);  % 形成矩阵，大小：Nrg*Naz

D_fn_ref_Vr = sqrt(1-lamda^2*fn_ref^2/(4*Vr^2));    % 参考频率fn_ref处的徙动因子，是常数。

K_src = 2*Vr^2*f0^3.*D_fn_Vr.^3./(c*R_ref*(fa.').^2);   % 列向量，使用R_ref处的值 
K_src_mtx = K_src*ones(1,Nrg);  % 形成矩阵
Km = Kr./(1-Kr./K_src_mtx);     % 矩阵，这是变换到距离多普勒域的距离调频率。
                                % 使用 R_ref 处的值

% 下面生成 变标方程 s_sc
s_sc = exp(1j*pi.*Km.*(D_fn_ref_Vr./D_fn_Vr_mtx-1).*(tau_mtx-2*R_ref./(c.*D_fn_Vr_mtx)).^2);

% 下面将距离多普勒域的信号与变标方程相乘，实现“补余RCMC”
S_RD_1 = S_RD.*s_sc;            % 相位相乘，实现“补余RCMC”

% 作图
figure;
imagesc(abs(S_RD));
title('原始数据变换到距离多普勒域，幅度');
figure;
imagesc(abs(S_RD_1));
title('距离多普勒域，补余RCMC后，幅度');

%% 
% --------------------------------------------------------------------
% 变换到二维频域，进行“距离压缩，SRC，一致RCMC”
% --------------------------------------------------------------------
S_2df_1 = fft(S_RD_1,NFFT_r,2);         % 进行距离向FFT，变换到二维频域。距离零频在两端

% 完成距离压缩，SRC，一致RCMC这三者相位补偿的滤波器为：
H1 = exp(1j*pi.*D_fn_Vr_mtx./(D_fn_ref_Vr.*Km).*fr_mtx.^2)...
    .*exp(1j*4*pi/c.*(1./D_fn_Vr_mtx-1/D_fn_ref_Vr).*R_ref.*fr_mtx);
% 上面的H1距离零频在中心
W_ref = ones(Naz,1)*(kaiser(Nrg,3).');	% 距离向，构建Kaiser窗，此为矩阵形式，距离零频在中心。
% H1 = W_ref.*H1;             % 加入距离平滑窗，以抑制旁瓣，距离零频在中心。
% 下面通过fftshift将H1的距离零频调整到两端
H1 = fftshift(H1,2);        % 左右半边互换，距离零频在两端。

S_2df_2 = S_2df_1.*H1;    	% 在二维频域，相位相乘，实现距离压缩，SRC，一致RCMC

S_RD_2 = ifft(S_2df_2,NFFT_r,2);    % 进行距离IFFT，回到距离多普勒域，完成所有距离处理。

% 作图
figure;
imagesc(abs(S_2df_1));
title('变换到二维频域');
figure;
imagesc(abs(S_2df_2));
title('相位相乘，实现距离压缩，SRC，一致RCMC后，二维频域');

figure;
imagesc(abs(S_RD_2));
title('完成距离压缩，SRC，一致RCMC后，距离多普勒域');

%%
% --------------------------------------------------------------------
% 距离多普勒域，完成“方位压缩”和“附加相位校正”
% --------------------------------------------------------------------
R0_RCMC = (c/2).*tau;   % 随距离线变化的R0，记为R0_RCMC，用来计算方位MF。
% 生成方位向匹配滤波器
Haz = exp(1j*4*pi.*(D_fn_Vr*R0_RCMC).*f0./c);       % 方位MF

% 附加相位校正项
H2 = exp(-1j*4*pi.*Km./(c^2).*(1-D_fn_Vr_mtx./D_fn_ref_Vr)...
    .*((1./D_fn_Vr)*R0_RCMC-R_ref./D_fn_Vr_mtx).^2); 	% 附加相位校正项

% 下面进行相位相乘，在距离多普勒域，同时完成方位MF和附加相位校正
S_RD_3 = S_RD_2.*Haz.*H2;           % 距离多普勒域，相位相乘

% 最后通过IFFT回到图像域，完成方未处理
s_image = ifft(S_RD_3,NFFT_a,1); 	% 完成成像过程，得到成像结果为：s_image

% 作图
figure;
imagesc(abs(S_RD_3));
title('距离多普勒域，进行了相位相乘后（方位MF和附加相位校正）');

figure;
imagesc(abs(s_image));
title('成像结果');

%%
% 下面通过调用函数，得到三个点目标各自的切片，并进行升采样
% 同时对点目标中心做距离向切片，方位向切片
% 计算出相应的指标：PSLR，ISLR，IRW
NN = 20;
% 分别得到每个点目标的切片放大；行切片、列切片；和相应的指标

% 目标1，点目标中心在 （ tg_1_x，tg_1_y ）
% =========================================================================
% 现在的点目标位置计算如下：
tg_1_x = rem( R0*tan(sita_r_c)/Vr*Fa , Naz );
if tg_1_x < Naz/2
    tg_1_x = tg_1_x + (Naz/2+1);
else
    tg_1_x = tg_1_x - (Naz/2+1);
end
tg_1_x = round(tg_1_x);    	% 四舍五入，得到整数值，作为点目标的方位中心坐标。
% 这里得到的 tg_1_x 即是点目标中心方位向的位置（坐标）。
% =========================================================================
% 下面计算目标1的距离向位置:
% 由于CSA的变标作用，从原来的压至零多普勒（R0），变为压至方位参考频率（fn_ref）处
% 的距离单元（即 R0/D_fn_ref_Vr ），因此对应的目标1的y轴位置如下，为 tg_1_y ：
tg_1_y = round( (Nrg/2+1) + 2*(R0/D_fn_ref_Vr-R0)/c*Fr );
target_1 = target_analysis_2( s_image(tg_1_x-NN:tg_1_x+NN,tg_1_y-NN:tg_1_y+NN),Fr,Fa,Vr);


% 目标2，点目标中心在 （tg_2_x，target_2_y）
tg_2_x = tg_1_x + delta_R1/Vr*Fa;
tg_2_y = tg_1_y;
% target_2 = target_analysis_2( s_image(tg_2_x-NN:tg_2_x+NN,tg_2_y-NN:tg_2_y+NN),Fr,Fa,Vr);


% 目标3，点目标中心在（tg_3_x，tg_3_y）
tg_3_x = tg_2_x + delta_R2*tan(sita_r_c)/Vr*Fa;
tg_3_x = fix(tg_3_x);
tg_3_y = tg_2_y + 2*(delta_R2/D_fn_ref_Vr)/c*Fr;
% target_3 = target_analysis_2( s_image(tg_3_x-NN:tg_3_x+NN,tg_3_y-NN:tg_3_y+NN),Fr,Fa,Vr);



