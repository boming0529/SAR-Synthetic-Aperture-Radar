%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            小斜視角（3.5°）
%                RDA
%      RDA in the Low Squint Case
%            點目標模擬
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update Nov 18 2023. 

%%
clear;
close all;
clc;
% 定義參數
% --------------------------------------------------------------------
% Table 6.1 : C-Band AirBorn SAR Parameters Used in the Simulations
% --------------------------------------------------------------------
R_nc = 20e3;                % 景中心斜距   (Slant range of scene center, Km)
Vr = 150;                   % 雷達有效速度 (Effective rader velocity, m/s)
Tr = 2.5e-6;                % 發射脈衝時寬 (Transmitted pulse duration, \mu s)
Kr = 20e12;                 % 距離向調频率 (Range FM rate, Hz/sec)
f0 = 5.3e9;                 % 雷達工作頻率 (Radar center frequency, GHz)
BW_dop = 80;                % 多普勒带寬   (Doppler bandwidth, Hz)
Fr = 60e6;                  % 距離採樣率   (Range sampling rate, MHz)
Fa = 100;                   % 方位採樣率   (Azimuth sampling rate or PRF, HZ)
Naz = 256;                  % 方位採樣數,即數據矩陣列數(Number of range lines(column))——這裡修改為1024。
Nrg = 320;                  % 距離採樣數,即數據矩阵行數(Sample per range line(row)）
sita_r_c = (3.5*pi)/180;	% 波束斜视角,3.5度,這裡轉為弧度(Beam squint angel, rad)

% --------------------------------------------------------------------
% 定義參數
% --------------------------------------------------------------------
c = 3e8;                    % 光速(Speed of light)
R0 = R_nc*cos(sita_r_c);	% 與R_nc相對應的最近斜距,根據(4.17)計算.(Range of closest approach)
Nr = Tr*Fr;                 % 線性調頻信號採樣點數,課本找不到
BW_range = Kr*Tr;           % 距離向帶寬.(表4.3).(Range bandwidth)
lamda = c/f0;               % 波長(Wavelength)
fnc = 2*Vr*sin(sita_r_c)/lamda;             %多普勒中心頻率,根據(4.33)計算.(Doppler centroid frequency)
La_real = 0.886*2*Vr*cos(sita_r_c)/BW_dop;  %方位向天線長度,根據(4.36)計算.(Antenna length)  
La = 0.886*R_nc*lamda/La_real;              %方位角波束覆蓋區.(表4.3).(Azimuth beam footprint)
a_sr = Fr / BW_range;       % 距離向過採樣率,根據(4.22)計算.(Range oversampling ratio) 
a_sa = Fa / BW_dop;         % 方位向過採樣率,尚未找到.(Azimuth oversampling factor)
Mamb = round(fnc/Fa);       % 多普勒模糊,根據(12.20)計算,與課本公式不同.(Doppler ambiguity)

NFFT_r = Nrg;               % 距離向FFT長度,即數據矩阵行數
NFFT_a = Naz;               % 方位向FFT長度,即數據矩陣列數

% 定義座標軸
% --------------------------------------------------------------------
% Figure 6.2 : Positions of three targets used in the simulation
% --------------------------------------------------------------------
% 模擬目標位置
% 以距離向作為x軸正方向
% 以方位向作為y軸正方向
% -------------------------------------------------------------------- 
delta_R0 = 0;       % 將目標A的波束中心穿越时刻，定義為方位向時間零點。
delta_R1 = 120; 	% 目標A和目標B的方位向距離差，120m
delta_R2 = 50;      % 目標B和目標C的距離向距离差，50m

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
nc_1 = (y1-x1*tan(sita_r_c))/Vr;    % 目標A的波束中心穿越时刻。
nc_2 = (y2-x2*tan(sita_r_c))/Vr;    % 目標B的波束中心穿越时刻。
nc_3 = (y3-x3*tan(sita_r_c))/Vr;    % 目標C的波束中心穿越时刻。
nc_target = [nc_1,nc_2,nc_3];       % 定義目標A,B,C穿越时刻矩陣。

% 定義tau和eta
% --------------------------------------------------------------------
% 距離(方位)向时间，频率相关定义
% --------------------------------------------------------------------
% 距离
tr = 2*x1/c + ( -Nrg/2 : (Nrg/2-1) )/Fr;                % 距離時間軸
fr = ( -NFFT_r/2 : NFFT_r/2-1 )*( Fr/NFFT_r );          % 距離頻率軸
% 方位
ta = ( -Naz/2: Naz/2-1 )/Fa;                            % 方位時間軸
fa = fnc + ( -NFFT_a/2 : NFFT_a/2-1 )*( Fa/NFFT_a );	% 方位頻率軸

% 生成距离（方位）时间（频率）矩阵
tau_mtx = ones(Naz,1)*tr;    % 距离时间轴矩阵，大小：Naz*Nrg
ta_mtx = ta.'*ones(1,Nrg);  % 方位时间轴矩阵，大小：Naz*Nrg

%% 
% --------------------------------------------------------------------
% 生成点目标原始数据
% --------------------------------------------------------------------
s_echo = zeros(Naz,Nrg);    % 用来存放生成的回波数据

A0 = 1;                     % 目标回波幅度，都设置为1.
for k = 1:3                 % 生成k个目标的原始回波数据
    R_n = sqrt( (x_range(k).*ones(Naz,Nrg)).^2 + (Vr.*ta_mtx-y_azimuth(k).*ones(Naz,Nrg)).^2 );% 目标k的瞬时斜距 (6.2)
    w_range = ((abs(tau_mtx-2.*R_n./c)) <= ((Tr/2).*ones(Naz,Nrg)));     % 距离向包络，即距离窗
    % =====================================================================    
    % 方位向包络，也就是 天线的双程方向图作用因子。
    %{
    % 方式1
    % sinc平方型函数，根据公式（4.31）计算
    beta_bw = 0.886*lamda/La_real;              %雷達3dB波束 (表4.3).(Azimuth beamwidth)
    sita = atan( Vr.*(ta_mtx-nc_target(k).*ones(Naz,Nrg))/x_range(k) );
    w_azimuth1 = (sinc(0.886.*sita./beta_bw)).^2;    
    % 用每个目标对应的 波束中心穿越时刻，而不是之前参数中的nc。
    %}
    %
    % 方式2
    % 利用合成孔径长度，直接构造矩形窗（其实这里只是限制数据范围，没有真正加窗）
    w_azimuth = (abs(ta - nc_target(k)) <= (La/2)/Vr);    % 行向量
    w_azimuth = w_azimuth.'*ones(1,Nrg);    % 生成Naz*Nrg的矩阵
    %}
    % =====================================================================     
    s_k = A0.*w_range.*w_azimuth.*exp(-(1j*4*pi*f0).*R_n./c).*exp((1j*pi*Kr).*(tau_mtx-2.*R_n./c).^2);  % (6.1)
    % s_k = A0.*w_range.*w_azimuth1.*w_azimuth.*exp(-(1j*4*pi*f0).*R_n./c).*exp((1j*pi*Kr).*(tau_mtx-2.*R_n./c).^2);
    % 上式就是生成的某一个点目标（目标k）的回波信号。

    % 经过几次循环，生成几个点目标的回波信号，相加即可。
    s_echo = s_echo + s_k;  % 所有目標點回波訊號之和   
end
% s_echo 就是我们需要的原始数据，点目标回波信号。

% 作圖
% 圖一： 低斜視角情況下的三點雷達原始模擬訊號 (課本圖 6.3)
figure;
subplot(2,2,1);
imagesc(real(s_echo));
title('(a) 實部 Real Part');
xlabel('Range Time Domain(Samples)');
ylabel('Azimuth Time Domain(Samples)');
text(280,-150,'圖一: 目標物 A 原始數據');     

subplot(2,2,2);
imagesc(imag(s_echo));
title('(b) 虚部 Imaginary Part');
xlabel('Range Time Domain(Samples)');
ylabel('Azimuth Time Domain(Samples)');

subplot(2,2,3);
imagesc(abs(s_echo));
title('(c) 強度 Magnitude');
xlabel('Range Time Domain(Samples)');
ylabel('Azimuth Time Domain(Samples)');

subplot(2,2,4);
imagesc(angle(s_echo));
title('(d) 相位 Phase');
xlabel('Range Time Domain(Samples)');
ylabel('Azimuth Time Domain(Samples)');
% colormap(gray);

%%
% --------------------------------------------------------------------
% 6.3.2 Range Compression
% --------------------------------------------------------------------
S_range = fft(s_echo,NFFT_r,2);     % 进行距离向傅里叶变换，零频在两端。

%
% 作圖
% 图2——距离频域，方位时域，频谱（未距离压缩）
figure;
subplot(1,2,1);
imagesc(real(S_range));
title('(a)實部');
xlabel('Range frequency Domain(Samples)');
ylabel('Azimuth time Domain(Samples)');
text(280,-60,'圖二: Range frequency Domain');       % 给图2进行文字说明
text(340,-10,'未壓縮');       

subplot(1,2,2);
imagesc(abs(S_range));
title('(b)幅度');
xlabel('Range frequency Domain(Samples)');
ylabel('Azimuth time Domain(Samples)');
%}

%　生成距离向匹配滤波器 (matched filter)
% ====================================================
% 采用方式2
% 时域复制脉冲，末端补零，fft，再取复共轭。
t_ref = ( -Nr/2 : (Nr/2-1) )/Fr;    % 用来生成距离MF的距离时间轴
t_ref_mtx = ones(Naz,1)*t_ref;      % 矩阵形式
w_ref = kaiser(Nr,2.5);             % 距离向，构建Kaiser窗，此为列向量。
w_ref = ones(Naz,1)*(w_ref.');      % 构成矩阵形式，每一行都相同的加窗。

% s_ref = exp((1j*pi*Kr).*((t_ref_mtx).^2)); % 复制（发射）脉冲，未加窗。
s_ref = w_ref.*exp((1j*pi*Kr).*((t_ref_mtx).^2)); % 复制（发射）脉冲，加了窗 (3.48)。

% region : (3.49)
s_ref = [s_ref,zeros(Naz,Nrg-Nr)];      % 对复制脉冲，后端补零。
S_ref = fft(s_ref,NFFT_r,2);            % 复制脉冲的距离傅里叶变换，零频在两端。
H_range = conj(S_ref);                  % 距离向匹配滤波器，零频在两端。
% endregion 

S_range_c = S_range.*H_range;           % 乘以匹配滤波器，零频在两端。    
s_rc = ifft(S_range_c,[],2);            % 完成距离压缩，回到二维时域。 (6.3)
% s_rc的长度为：Naz*Nrg。未去除弃置区。

% 对s_rc进行去除弃置区的操作
% 弃置区长度为：2*（Nr-1）
% 我们截取的长度：（Nrg-Nr+1），记为 N_rg。
N_rg = Nrg-Nr+1;                        % 完全卷积的长度
s_rc_c = zeros(Naz,N_rg);               % 用来存放去除弃置区后的数据
s_rc_c = s_rc(:,1:N_rg);                % 取前 N_rg列。
% ====================================================

%
% 作圖
% 图3——距离频域，方位时域，频谱（已距离压缩）
figure;
subplot(1,2,1);
imagesc(real(S_range_c));
title('(a)實部');
xlabel('Range frequency Domain(Samples)');
ylabel('Azimuth time Domain(Samples)');
text(280,-60,'圖三: Range frequency Domain');       % 给图3进行文字说明
text(340,-10,'已壓縮');       

subplot(1,2,2);
imagesc(abs(S_range_c));
title('(b)幅度');
xlabel('Range frequency Domain(Samples)');
ylabel('Azimuth time Domain(Samples)');
%}
%

% 作圖
% 圖四 :Range Compressed Data (課本圖 6.4)
figure;
subplot(1,2,1);
imagesc(real(s_rc_c));  %　这及其以下，都直接使用去除弃置区后的结果
title('(a)實部');
xlabel('Range Time Domain(Samples)');
ylabel('Azimuth Time Domain(Samples)');
text(130, -15,'圖四 : Range Compressed Data');       % 完成壓縮

subplot(1,2,2);
imagesc(abs(s_rc_c));
title('(b)幅度');
xlabel('Range Time Domain(Samples)');
ylabel('Azimuth Time Domain(Samples)');
%}

%%
% --------------------------------------------------------------------
% 6.3.4 RCMC
% 变换到距离多普勒域，进行距离徙动校正 
% --------------------------------------------------------------------
s_rc_c = s_rc_c.*exp(-1j*2*pi*fnc.*(ta.'*ones(1,N_rg)));    % 数据搬移
S_rd = fft(s_rc_c,NFFT_a,1);            % 方位向傅里叶变换，到距离多普勒域 (Azimuth FFT)
% ====================================================================
% 设置方位向频率轴
fa = fnc + fftshift(-NFFT_a/2:NFFT_a/2-1)/NFFT_a*Fa;    % 方位频率轴如此设置。
% =====================================================================
% 下面这个是改进的，每一个最近斜距（R0）都随着距离门的不同而改变。
tr_RCMC = 2*x1/c + ( -N_rg/2 : (N_rg/2-1) )/Fr;   % 在新的距离线长度下的时间轴。

R0_RCMC = (c/2).*tr_RCMC*cos(sita_r_c);   % 随距离线变化的R0，记为R0_RCMC，用来计算RCM和Ka。
delta_Rrd_fn = lamda^2.*((fa.').^2)*(R0_RCMC)/(8*Vr^2);

num_range = c/(2*Fr);   % 一个距离采样单元，对应的长度
delta_Rrd_fn_num = delta_Rrd_fn./num_range; % 每一个方位向频率，其RCM对应的距离采样单元数

R = 8;  % sinc插值核长度
S_rd_rcmc = zeros(NFFT_a,N_rg); % 用来存放RCMC后的值
for p = 1 : NFFT_a
    for q = 1 : N_rg   % 此时距离向的长度是 (Nrg-Nr+1)=N_rg        
        delta_Rrd_fn_p = delta_Rrd_fn_num(p,q);
        Rrd_fn_p = q + delta_Rrd_fn_p;
        Rrd_fn_p_zheng = ceil(Rrd_fn_p);        % ceil，向上取整。
        ii = ( Rrd_fn_p-(Rrd_fn_p_zheng-R/2):-1:Rrd_fn_p-(Rrd_fn_p_zheng+R/2-1)  );        
        rcmc_sinc = sinc(ii);
        rcmc_sinc = rcmc_sinc/sum(rcmc_sinc);   % 插值核的归一化
        % ii 是sinc插值过程的变量; (2.7)
        % g(x)=sum(h(ii)*g_d(x-ii)) = sum(h(ii)*g_d(ll));
               
        % 由于S_rd只有整数点取值，且范围有限。因此插值中要考虑它的取值溢出边界问题。
        % 这里我采取循环移位的思想，用来解决取值溢出问题。
        if (Rrd_fn_p_zheng-R/2) > N_rg    % 全右溢
            ll = (Rrd_fn_p_zheng-R/2-N_rg:1:Rrd_fn_p_zheng+R/2-1-N_rg);
        else
            if (Rrd_fn_p_zheng+R/2-1) > N_rg    % 部分右溢
                ll_1 = (Rrd_fn_p_zheng-R/2:1:N_rg);
                ll_2 = (1:1:Rrd_fn_p_zheng+R/2-1-N_rg);
                ll = [ll_1,ll_2];
            else
                if (Rrd_fn_p_zheng+R/2-1) < 1    % 全左溢（不可能发生，但还是要考虑）
                    ll = (Rrd_fn_p_zheng-R/2+N_rg:1:Rrd_fn_p_zheng+R/2-1+N_rg);
                else
                    if (Rrd_fn_p_zheng-R/2) < 1       % 部分左溢
                        ll_1 = (Rrd_fn_p_zheng-R/2+N_rg:1:N_rg);
                        ll_2 = (1:1:Rrd_fn_p_zheng+R/2-1);
                        ll = [ll_1,ll_2];
                    else
                        ll = (Rrd_fn_p_zheng-R/2:1:Rrd_fn_p_zheng+R/2-1);
                    end                    
                end
            end
        end   
        rcmc_S_rd = S_rd(p,ll);
        S_rd_rcmc(p,q) = sum( rcmc_sinc.*rcmc_S_rd );
    end
end
% S_rd_rcmc 就是RCMC后的距离多普勒域频谱。

% 作圖
% 圖5 : Range Dopper Domain（未RCMC）備註：課本沒有這張圖
figure;
subplot(1,2,1);
imagesc(real(S_rd));
title('(a)實部');
xlabel('Range time Domain(Samples)');
ylabel('Azimuth frequency Domain(Samples)');
text(130,-10,'圖五 : Range Dopper Domain');   % 未RCMC

subplot(1,2,2);
imagesc(abs(S_rd));
title('(b)幅度');
xlabel('Range time Domain(Samples)');
ylabel('Azimuth frequency Domain(Samples)');

% 作圖
% 圖六： RCMC 後的模擬數據 (課本圖 6.9)
figure;
subplot(1,2,1);
imagesc(real(S_rd_rcmc));
title('(a)實部');
xlabel('Range time Domain(Samples)');
ylabel('Azimuth frequency Domain(Samples)');
text(130,-10,'圖六 : RCMC 後的模擬數據');       % 已RCMC  

subplot(1,2,2);
imagesc(abs(S_rd_rcmc));
title('(b)幅度');
xlabel('Range time Domain(Samples)');
ylabel('Azimuth frequency Domain(Samples)');
%}

%%
% --------------------------------------------------------------------
% 6.3.6 Azimuth Compression
% --------------------------------------------------------------------
fa_azimuth_MF = fa;         % 方位频率轴，采用和RCMC中所用的频率轴相同。
Ka = 2*Vr^2*(cos(sita_r_c))^3./(lamda.* R0_RCMC);  	% 方位向调频率，是随最近斜距R0变化的。
Ka_1 = 1./Ka;                                       % 为了方便计算，先取倒数。
Haz = exp( -1j*pi.*(((fa_azimuth_MF).').^2*Ka_1) );	% 方位向匹配滤波器
% 这里要注意，生成的MF的零频既不是在两端，也不是在中心的。
% 考虑下频率轴是什么样的，间断点在哪里。注意fa的构成。
% 这里的频率轴和距离多普勒域的方位频谱是对应的。

S_rd_c = S_rd_rcmc.*Haz;            % 乘以匹配滤波器
s_ac = ifft(S_rd_c,[],1);       	% 完成方位压缩，变到图像域。结束。

% 作图
% 图7——成像结果
figure;
imagesc(abs(s_ac));
title('点目标成像');
xlabel('距离时域（采样点）');
ylabel('方位时域（采样点）');     

%%
% 下面通过调用函数，得到三个点目标各自的切片，并进行升采样
% 同时对点目标中心做距离向切片，方位向切片
% 计算出相应的指标：PSLR，ISLR，IRW

NN = 20;
% 分别得到每个点目标的切片放大；行切片、列切片；和相应的指标
% 目标1，点目标中心在 （ tg_1_x，tg_1_y ）
tg_1_x = 96;
tg_1_y = round(N_rg/2);
% target_1 = target_analysis_2( s_ac(tg_1_x-NN:tg_1_x+NN,tg_1_y-NN:tg_1_y+NN),Fr,Fa,Vr);

% 目标2，点目标中心在 （tg_2_x，target_2_y）
tg_2_x = tg_1_x + delta_R1/Vr*Fa;
tg_2_y = tg_1_y;
% target_2 = target_analysis_2( s_ac(tg_2_x-NN:tg_2_x+NN,tg_2_y-NN:tg_2_y+NN),Fr,Fa,Vr);

% 目标3，点目标中心在（tg_3_x，tg_3_y）
tg_3_x = tg_2_x + delta_R2*tan(sita_r_c)/Vr*Fa;
tg_3_x = fix(tg_3_x);
tg_3_y = tg_2_y + 2*delta_R2/c*Fr;
% target_3 = target_analysis_2( s_ac(tg_3_x-NN:tg_3_x+NN,tg_3_y-NN:tg_3_y+NN),Fr,Fa,Vr);


