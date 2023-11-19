%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            小斜視角（3.5°）
%                RDA
%      RDA in the Low Squint Case
%            點目標模擬
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update Nov 18 2023. 

% --------------------------------------------------------------------------------------------------------
% 清除workspace和figure以及Command Window
% --------------------------------------------------------------------------------------------------------
clear;
close all;
clc;

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
Fa = 100;                   % 方位採樣率   (Azimuth sampling rate or PRF, HZ)
Naz = 256;                  % 方位採樣數,即數據矩陣列數(Number of range lines(column))——這裡修改為1024。(◕‿◕)?
Nrg = 320;                  % 距離採樣數,即數據矩陣行數(Sample per range line(row)）
sita_r_c = (3.5*pi)/180;	% 波束斜视角,3.5度,這裡轉為弧度(Beam squint angel, rad)

% --------------------------------------------------------------------------------------------------------
% 定義參數
% --------------------------------------------------------------------------------------------------------
c = 3e8;                    % 光速(Speed of light)
R0 = R_nc*cos(sita_r_c);	% 與R_nc相對應的最近斜距,根據(4.17)計算.(Range of closest approach)
Nr = Tr*Fr;                 % 線性調頻信號採樣點數,課本找不到 (◕‿◕)?
BW_range = Kr*Tr;           % 距離向帶寬.(表4.3).(Range bandwidth)
lamda = c/f0;               % 波長(Wavelength)
fnc = 2*Vr*sin(sita_r_c)/lamda;             %多普勒中心頻率,根據(4.33)計算.(Doppler centroid frequency)
La_real = 0.886*2*Vr*cos(sita_r_c)/BW_dop;  %方位向天線長度,根據(4.36)計算.(Antenna length)  
La = 0.886*R_nc*lamda/La_real;              %方位角波束覆蓋區.(表4.3).(Azimuth beam footprint)
a_sr = Fr / BW_range;       % 距離向過採樣率,根據(4.22)計算.(Range oversampling ratio) 
a_sa = Fa / BW_dop;         % 方位向過採樣率,尚未找到.(Azimuth oversampling factor)
Mamb = round(fnc/Fa);       % 多普勒模糊,根據(12.20)計算,與課本公式不同.(Doppler ambiguity)

NFFT_a = Naz;               % 方位向FFT長度,即數據矩陣列數

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
nc_1 = (y1-x1*tan(sita_r_c))/Vr;    % 目標A的波束中心穿越时刻
nc_2 = (y2-x2*tan(sita_r_c))/Vr;    % 目標B的波束中心穿越时刻
nc_3 = (y3-x3*tan(sita_r_c))/Vr;    % 目標C的波束中心穿越时刻
nc_target = [nc_1,nc_2,nc_3];       % 定義目標A,B,C穿越时刻矩陣

% --------------------------------------------------------------------------------------------------------
% 定義tau和eta矩陣
% --------------------------------------------------------------------------------------------------------
% 距離列矩陣
tau = 2*x1/c + ( -Nrg/2 : (Nrg/2-1) )/Fr;  % 距離時間軸,課本找不到(◕‿◕)
% 生成距離矩陣
tau_mtx = ones(Naz,1)*tau;                 % 距離時間軸矩陣,矩陣大小：Naz*Nrg

% 方位列矩陣
eta = ( -Naz/2: Naz/2-1 )/Fa;              % 方位時間軸,課本找不到(◕‿◕)
% 生成方位矩陣
eta_mtx = eta.'*ones(1,Nrg);               % 方位時間轴矩陣,矩陣大小：Naz*Nrg

% --------------------------------------------------------------------------------------------------------
% 6.3.1 Raw data
% --------------------------------------------------------------------------------------------------------
s_echo = zeros(Naz,Nrg);    % s_echo用来存放生成的回波(Raw Data)結果

A0 = 1;                     % 目標回波幅度,都設置為1
for k = 1:3                 % 生成k個目標的原始回波數據
    R_n = sqrt( (x_range(k).*ones(Naz,Nrg)).^2 + (Vr.*eta_mtx-y_azimuth(k).*ones(Naz,Nrg)).^2 );% 目標k的瞬時斜距(6.2)
    w_range = ((abs(tau_mtx-2.*R_n./c)) <= ((Tr/2).*ones(Naz,Nrg)));                            % 距離向包络,即距離窗(Range evelope)
    % ===================================================================================================   
    % 利用合成孔径长度,直接构造矩形窗(其实这里只是限制数据范围,没有真正加窗）(◕‿◕)?

    w_azimuth = (abs(eta - nc_target(k)) <= (La/2)/Vr);    % 方位向包络,行向量(Azimuth evelope) 
    w_azimuth = w_azimuth.'*ones(1,Nrg);                   % 生成Naz*Nrg的矩陣
    % ===================================================================================================    
    % 下式就是生成某一個點目標(目標k)的回波信號

    s_k = A0.*w_range.*w_azimuth.*exp(-(1j*4*pi*f0).*R_n./c).*exp((1j*pi*Kr).*(tau_mtx-2.*R_n./c).^2);  % (6.1)
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

% --------------------------------------------------------------------------------------------------------
% 6.3.2 Range Compression
% --------------------------------------------------------------------------------------------------------
% Range FFT 
% --------------------------------------------------------------------------------------------------------
% Make Figure 2, book doesn't have this figure  
% --------------------------------------------------------------------------------------------------------
NFFT_r = Nrg;                       % 距離向FFT長度,即數據矩陣行數
S_range = fft(s_echo,NFFT_r,2);     % 進行距離向FFT轉換,零頻在兩端

figure('Name','Range FFT', 'NumberTitle','on'); 

subplot(1,2,1);
imagesc(real(S_range));
title('(a)實部');
xlabel('Range frequency domain(Samples)');
ylabel('Azimuth time domain(Samples)');
   
subplot(1,2,2);
imagesc(abs(S_range));
title('(b)強度');
xlabel('Range frequency domain(Samples)');
ylabel('Azimuth time domain(Samples)');

% --------------------------------------------------------------------------------------------------------
%　Range matched filter
% --------------------------------------------------------------------------------------------------------
% Matched filter(MF)採用課本P55方式2:Range複製脈衝,末端補零,離散FFT,再取複共軛。

% 運算(3.48)
t_ref = ( -Nr/2 : (Nr/2-1) )/Fr;    % 用來生成MF的距離時間軸,課本找不到       (◕‿◕)?
t_ref_mtx = ones(Naz,1)*t_ref;      % 構建MF的距離時間軸矩陣,矩陣大小：Naz*Nr (◕‿◕)?
pulse_replia = exp((1j*pi*Kr).*((t_ref_mtx).^2));  % 複製脈衝訊號,根據(3.48)計算.(pulse replia)
w_ref = kaiser(Nr,2.5);             % 構建Kaiser窗,此为列向量,課本P56,w(t,T).(Tapering window)
w_ref = ones(Naz,1)*(w_ref.');      % 構建矩陣形式,每一行都相同的加窗,矩陣大小：Naz*Nr

h_prime_t = w_ref.*pulse_replia;    % 複製脈衝訊號,構建Kaiser窗(3.48)

% 運算(3.49)
h_prime_t = [h_prime_t,zeros(Naz,Nrg-Nr)]; % 對複製脈衝訊號,後端補零,稱為複製脈衝補零訊號
H2_f = fft(h_prime_t,NFFT_r,2);            % 對複製脈衝補零訊號取距離FFT,零频在兩端
H_range = conj(H2_f);                      % 取距離FFT後的複數共軛

S_range_c = S_range.*H_range;              % 在距離向頻率域將訊號壓縮

% --------------------------------------------------------------------------------------------------------
% Make Figure 3 : book doesn't have this figure
% --------------------------------------------------------------------------------------------------------
figure('Name','Range FFT X Matched Filter', 'NumberTitle','on'); 
subplot(1,2,1);
imagesc(real(S_range_c));
title('(a)實部');
xlabel('Range frequency domain(Samples)');
ylabel('Azimuth time domain(Samples)');

subplot(1,2,2);
imagesc(abs(S_range_c));
title('(b)強度');
xlabel('Range frequency domain(Samples)');
ylabel('Azimuth time domain(Samples)');

% --------------------------------------------------------------------------------------------------------
% Get Range compresses date 
% --------------------------------------------------------------------------------------------------------
% 運算(6.3)
s_rc = ifft(S_range_c,[],2);            % 完成距離向壓縮後,IFFT距離向回到時域,矩陣大小Naz*Nrg (6.3)
% 取得Range compresses date 
N_rg = Nrg-Nr+1;                        % 完全卷積的长度=截取的长度(Nrg-Nr+1)標記為,N_rg,矩陣大小Naz*N_rg
s_rc_c = s_rc(:,1:N_rg);                % :表示row全拿,第一行拿到N_rg行,課本沒有,此步驟才可完成Range compresses date 

% --------------------------------------------------------------------------------------------------------
%  Make Figure 4 : Range compresses date 
% --------------------------------------------------------------------------------------------------------
figure('Name','Range Compressed Data', 'NumberTitle','on');
subplot(1,2,1);
imagesc(real(s_rc_c));  
title('(a)實部');
xlabel('Range Time domain(Samples)');
ylabel('Azimuth Time domain(Samples)');

subplot(1,2,2);
imagesc(abs(s_rc_c));
title('(b)強度');
xlabel('Range Time domain(Samples)');
ylabel('Azimuth Time domain(Samples)');

% --------------------------------------------------------------------------------------------------------
% 6.3.3 Azimuth Fourier Transform
% --------------------------------------------------------------------------------------------------------
% 變換到距離多普勒域 
% --------------------------------------------------------------------------------------------------------
s_rc_c = s_rc_c.*exp(-1j*2*pi*fnc.*(eta.'*ones(1,N_rg)));    % 數據搬移(◕‿◕)?
S_rd = fft(s_rc_c,NFFT_a,1);            % 方位向FFT變化,到距離多普勒域(Azimuth FFT)

% --------------------------------------------------------------------------------------------------------
%  Make Figure 5 : Azimuth FFT of the low squint example
% --------------------------------------------------------------------------------------------------------
figure('Name','Azimuth FFT of the low squint example');
subplot(1,2,1);
imagesc(real(S_rd));
title('(a)實部');
xlabel('Range time domain(Samples)');
ylabel('Azimuth frequency domain(Samples)'); 

subplot(1,2,2);
imagesc(abs(S_rd));
title('(b)強度');
xlabel('Range time domain(Samples)');
ylabel('Azimuth frequency domain(Samples)');

% --------------------------------------------------------------------------------------------------------
% 6.3.4 RCMC
% --------------------------------------------------------------------------------------------------------
% 進行距離徙動校正 
% --------------------------------------------------------------------------------------------------------

% 運算(6.11)
% 設置方位向頻率軸
fn = fnc + fftshift(-NFFT_a/2:NFFT_a/2-1)/NFFT_a*Fa;        % 方位频率轴

% 下面这个是改进的，每一个最近斜距（R0）都随着距离门的不同而改变。
tr_RCMC = 2*x1/c + ( -N_rg/2 : (N_rg/2-1) )/Fr;            % 在新的距离线长度下的时间轴。

R0_RCMC = (c/2).*tr_RCMC*cos(sita_r_c);                    % 随距离线变化的R0，记为R0_RCMC，用来计算RCM和Ka。
delta_Rrd_fn = lamda^2.*((fn.').^2)*(R0_RCMC)/(8*Vr^2);

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
% 圖六： RCMC 後的模擬數據 (課本圖 6.9)
figure;
subplot(1,2,1);
imagesc(real(S_rd_rcmc));
title('(a)實部');
xlabel('Range time domain(Samples)');
ylabel('Azimuth frequency domain(Samples)');
text(130,-10,'圖六 : RCMC 後的模擬數據');       % 已RCMC  

subplot(1,2,2);
imagesc(abs(S_rd_rcmc));
title('(b)強度');
xlabel('Range time domain(Samples)');
ylabel('Azimuth frequency domain(Samples)');
%}

%%
% --------------------------------------------------------------------
% 6.3.6 Azimuth Compression
% --------------------------------------------------------------------
fa_azimuth_MF = fn;         % 方位频率轴，采用和RCMC中所用的频率轴相同。
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


