%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            小斜視角（3.5°）
%                RDA
%      RDA in the Low Squint Case
%            點目標模擬
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update Nov 18 2023. 

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
Fa = 200;                   % 方位採樣率   (Azimuth sampling rate or PRF, HZ)——這裡修改為200.(◕‿◕)?
Naz = 1024;                  % 方位採樣數,即數據矩陣列數(Number of range lines(row))——這裡修改為1024.(◕‿◕)?
Nrg = 320;                  % 距離採樣數,即數據矩陣行數(Sample per range line(column)）
sita_r_c = (3.5*pi)/180;	% 波束斜视角,3.5度,這裡轉為弧度(Beam squint angel, rad)

% --------------------------------------------------------------------------------------------------------
% 定義參數
% --------------------------------------------------------------------------------------------------------
c = 3e8;                    % 光速(Speed of light)
R0 = R_nc*cos(sita_r_c);	% 與R_nc相對應的最近斜距,根據(4.17)計算.(Range of closest approach)
Nr = Tr*Fr;                 % 線性調頻信號採樣點數,課本找不到(◕‿◕)?
BW_range = Kr*Tr;           % 距離向帶寬.(表4.3).(Range bandwidth)
lamda = c/f0;               % 波長(Wavelength)
fnc = 2*Vr*sin(sita_r_c)/lamda;             %多普勒中心頻率,根據(4.33)計算.(Doppler centroid frequency)
La_real = 0.886*2*Vr*cos(sita_r_c)/BW_dop;  %方位向天線長度,根據(4.36)計算.(Antenna length)  
La = 0.886*R_nc*lamda/La_real;              %方位角波束覆蓋區.(表4.3).(Azimuth beam footprint)
a_sr = Fr / BW_range;       % 距離向過採樣率,根據(4.22)計算.(Range oversampling ratio) 
a_sa = Fa / BW_dop;         % 方位向過採樣率,尚未找到.(Azimuth oversampling factor)
Mamb = round(fnc/Fa);       % 多普勒模糊,根據(12.20)計算,與課本公式不同(◕‿◕)?.(Doppler ambiguity)

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

%%
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
    % 利用合成孔徑長度,直接構造矩形窗(其實這里只是限制數據範圍,沒有真正加窗）(◕‿◕)?

    w_azimuth = (abs(eta - nc_target(k)) <= (La/2)/Vr);    % 方位向包络(Azimuth evelope) 
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

%%
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
w_ref = kaiser(Nr,2.5);             % 構建Kaiser窗,此為列向量,課本P56,w(t,T).(Tapering window)
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
s_rc_c = s_rc(:,1:N_rg);                % :表示row全拿,第一行拿到N_rg行,課本找不到(◕‿◕)?,此步驟才可完成Range compresses date 

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
s_rc_c = s_rc_c.*exp(-1j*2*pi*fnc.*(eta.'*ones(1,N_rg)));    % 數據搬移,課本找不到(◕‿◕)?
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
% 定義RCMC tau和RCMC eta矩陣
% --------------------------------------------------------------------------------------------------------
tau_RCMC = 2*x1/c + ( -N_rg/2 : (N_rg/2-1) )/Fr;            % 定義壓縮後的距離時間軸,類似程式89行
f_eta = fnc + fftshift(-NFFT_a/2:NFFT_a/2-1)/NFFT_a*Fa;     % 定義壓縮後的方位频率轴(◕‿◕)?

% --------------------------------------------------------------------------------------------------------

% 運算R0
R0_RCMC = (c/2).*tau_RCMC*cos(sita_r_c);                    % 由三角函數計算出R0

% 運算(6.11)
delta_Rrd_fn = lamda^2.*((f_eta.').^2)*(R0_RCMC)/(8*Vr^2);  % 校正的距離徒動,表示距離多普勒域中心多遠 
number_range = c/(2*Fr);                       % 計算壓縮後的距離採樣數,P149(sample spacing)
delta_Rrd_fn_num = delta_Rrd_fn./number_range; % 校正的距離徒動/計算壓縮後的距離採樣數,逼近距離多普勒域中心

% 運算插值(◕‿◕)?(◕‿◕)?(◕‿◕)?(◕‿◕)?
% 利用h(x)=sinc(x)跟gd(i)來做RMCM校正,找到多普勒中心
R = 8;                          % sinc插值該長度
S_rd_rcmc = zeros(NFFT_a,N_rg); % 用來存放RCMC後的值
for p = 1 : NFFT_a
    for q = 1 : N_rg            % 此时距離向的長度是 (Nrg-Nr+1)=N_rg        
        delta_Rrd_fn_p = delta_Rrd_fn_num(p,q);
        Rrd_fn_p = q + delta_Rrd_fn_p;
        Rrd_fn_p_zheng = ceil(Rrd_fn_p);        % ceil,向上取整
        ii = ( Rrd_fn_p-(Rrd_fn_p_zheng-R/2):-1:Rrd_fn_p-(Rrd_fn_p_zheng+R/2-1)  );        
        rcmc_sinc = sinc(ii);                   % ii 是sinc插值過程的變量; (2.7)
        rcmc_sinc = rcmc_sinc/sum(rcmc_sinc);   % 插值捲積核的歸一化(interpolation kernal)

% 處理校正後矩陣大小問題

        % 由於S_rd只有整數點取值，且範圍有限。因此插值中要考慮它的取值溢出邊界問題
        % 采取循環移位的思想，用來解決取值溢出問題
        if (Rrd_fn_p_zheng-R/2) > N_rg          % 全右溢
            ll = (Rrd_fn_p_zheng-R/2-N_rg:1:Rrd_fn_p_zheng+R/2-1-N_rg);
        else
            if (Rrd_fn_p_zheng+R/2-1) > N_rg    % 部分右溢
                ll_1 = (Rrd_fn_p_zheng-R/2:1:N_rg);
                ll_2 = (1:1:Rrd_fn_p_zheng+R/2-1-N_rg);
                ll = [ll_1,ll_2];
            else
                if (Rrd_fn_p_zheng+R/2-1) < 1    % 全左溢（不可能發生，但還是要考慮）
                    ll = (Rrd_fn_p_zheng-R/2+N_rg:1:Rrd_fn_p_zheng+R/2-1+N_rg);
                else
                    if (Rrd_fn_p_zheng-R/2) < 1  % 部分左溢
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
 % 運算(2.57)       
        S_rd_rcmc(p,q) = sum( rcmc_sinc.*rcmc_S_rd );
    end
end
% S_rd_rcmc 就是RCMC後的距離多普勒域頻譜

% --------------------------------------------------------------------------------------------------------
%  Make Figure 6 : Simulated data after RCMC
% --------------------------------------------------------------------------------------------------------
figure('Name', 'Simulated data after RCMC', 'NumberTitle', 'on')
subplot(1,2,1);
imagesc(real(S_rd_rcmc));
title('(a)實部');
xlabel('Range time domain(Samples)');
ylabel('Azimuth frequency domain(Samples)');   

subplot(1,2,2);
imagesc(abs(S_rd_rcmc));
title('(b)強度');
xlabel('Range time domain(Samples)');
ylabel('Azimuth frequency domain(Samples)');
% --------------------------------------------------------------------------------------------------------
% 6.3.6 Azimuth Compression
% --------------------------------------------------------------------------------------------------------
% Matched filter(MF)採用課本P55方式3:直接在頻率域產生匹配濾波器

fa_azimuth_MF = f_eta;                              % 方位頻率軸,採用和RCMC中所用的頻率軸相同
Ka = 2*Vr^2*(cos(sita_r_c))^3./(lamda.* R0_RCMC);  	% 計算Ka(4.38)
Ka_1 = 1./Ka;                                       % 為了方便計算,先取倒數
Haz = exp( -1j*pi.*(((fa_azimuth_MF).').^2*Ka_1) );	% 方位向匹配濾波器(6.16)
% 這裡要注意,生成的MF的零頻既不是在兩端,也不是在中心
% 考慮頻率軸是什麽樣的,間斷點在哪裡,注意feta的構成
% 這裡的頻率軸和距離多普勒域的方位頻譜是對應的

S_rd_c = S_rd_rcmc.*Haz;            % 乘以匹配濾波器(6.17)
s_ac = ifft(S_rd_c,[],1);       	% 完成方位向壓縮後,IFFT方位向回到時域(6.18)

% --------------------------------------------------------------------------------------------------------
%  Make Figure 7 : AzimuthCompressed Signal
% --------------------------------------------------------------------------------------------------------
figure('Name','Azimuth Compressed Signal', 'NumberTitle', 'on');
imagesc(abs(s_ac));
title('目標A,B,C成像');
xlabel('Range Time domain(Samples)');
ylabel('Azimuth Time domain(Samples)');    

% --------------------------------------------------------------------------------------------------------
%  Make Figure 8 : images power
% --------------------------------------------------------------------------------------------------------
figure('Name','Image Power', 'NumberTitle', 'on');
[row, col] = size(s_ac);
[X, Y] = meshgrid(1:col, 1:row);
mesh(X, Y, abs(s_ac));
title('point target A,B,C image power');
xlabel('Range Time domain(Samples)');
ylabel('Azimuth Time domain(Samples)');  
% --------------------------------------------------------------------------------------------------------


%%
% 下面通过调用函数，得到三个点目标各自的切片，并进行升采样
% 同时对点目标中心做距离向切片，方位向切片
% 计算出相应的指标：PSLR，ISLR，IRW

NN = 20;
% 分别得到每个点目标的切片放大；行切片、列切片；和相应的指标
% 目标1，点目标中心在 （ target_A_azimuth_samples, target_A_range_samples ）
target_A_azimuth_samples = 96;
target_A_range_samples = round(N_rg/2);
% target C 
% NN = 7;
% target_A_azimuth_samples = 130;
% target_A_range_samples = 106;

% --------------------------------------------------------------------------------------------------------
%  Make Figure 9 : point target C images power
% --------------------------------------------------------------------------------------------------------
figure('Name','Image Power', 'NumberTitle', 'on');
point_target_a_of_s_ac = s_ac(target_A_azimuth_samples-NN:target_A_azimuth_samples+NN,target_A_range_samples-NN:target_A_range_samples+NN); 
[row_a, col_a] = size(point_target_a_of_s_ac);
[X_a, Y_a] = meshgrid(1:col_a, 1:row_a);
mesh(target_A_range_samples-NN+X_a-1, target_A_azimuth_samples-NN+Y_a-1, abs(point_target_a_of_s_ac));
title('point target A image power');
xlabel('Range Time domain(Samples)');
ylabel('Azimuth Time domain(Samples)');  
% --------------------------------------------------------------------------------------------------------

% target_1 = target_analysis_2( s_ac(target_A_azimuth_samples-NN:target_A_azimuth_samples+NN,target_A_range_samples-NN:target_A_range_samples+NN),Fr,Fa,Vr);



% 目标B，点目标中心在 （ target_B_azimuth_samples, target_B_range_samples ）
target_B_azimuth_samples = target_A_azimuth_samples + delta_R1/Vr*Fa;
target_B_range_samples = target_A_range_samples;
% target_2 = target_analysis_2( s_ac(target_B_azimuth_samples-NN:target_B_azimuth_samples+NN,target_B_range_samples-NN:target_B_range_samples+NN),Fr,Fa,Vr);

% 目标C，点目标中心在 （ target_C_azimuth_samples, target_C_range_samples ）
target_C_azimuth_samples = target_B_azimuth_samples + delta_R2*tan(sita_r_c)/Vr*Fa;
target_C_azimuth_samples = fix(target_C_azimuth_samples);
target_C_range_samples = target_B_range_samples + 2*delta_R2/c*Fr;

% --------------------------------------------------------------------------------------------------------
%  Make Figure 10 : Expanded Target C
% --------------------------------------------------------------------------------------------------------
chip_size = 8
figure('Name','Expanded Target C', 'NumberTitle', 'on');
expanded_target_c_of_s_ac = s_ac(target_C_azimuth_samples-chip_size:target_C_azimuth_samples+chip_size,target_C_range_samples-chip_size:target_C_range_samples+chip_size); 
imagesc(abs(expanded_target_c_of_s_ac));
title('expanded target C image power');
xlabel('Range Time domain(Samples)');
ylabel('Azimuth Time domain(Samples)');  
% --------------------------------------------------------------------------------------------------------


% --------------------------------------------------------------------------------------------------------
%  Make Figure 10 : Expanded Target C (3D)
% --------------------------------------------------------------------------------------------------------
chip_size = 8
figure('Name','Expanded Target C (3D)', 'NumberTitle', 'on');
expanded_target_c_of_s_ac = s_ac(target_C_azimuth_samples-chip_size:target_C_azimuth_samples+chip_size,target_C_range_samples-chip_size:target_C_range_samples+chip_size); 
[row_expanded_c, col_expanded_c] = size(expanded_target_c_of_s_ac);
[X_expanded_c, Y_expanded_c] = meshgrid(1:col_expanded_c, 1:row_expanded_c);
mesh(target_C_range_samples-chip_size+ X_expanded_c -1, target_C_azimuth_samples-chip_size+Y_expanded_c -1, abs(expanded_target_c_of_s_ac));
title('expanded target C image power');
xlabel('Range Time domain(Samples)');
ylabel('Azimuth Time domain(Samples)');  
% --------------------------------------------------------------------------------------------------------

target_3 = target_analysis_2( s_ac(target_C_azimuth_samples-NN:target_C_azimuth_samples+NN,target_C_range_samples-NN:target_C_range_samples+NN),Fr,Fa,Vr);

% --------------------------------------------------------------------------------------------------------
%  Make Figure 11 : point target C images power
% --------------------------------------------------------------------------------------------------------
figure('Name','Image Power', 'NumberTitle', 'on');
point_target_c_of_s_ac = s_ac(target_C_azimuth_samples-NN:target_C_azimuth_samples+NN,target_C_range_samples-NN:target_C_range_samples+NN); 
[row_c, col_c] = size(point_target_c_of_s_ac);
[X_c, Y_c] = meshgrid(1:col_c, 1:row_c);
mesh(target_C_range_samples-NN+ X_c -1, target_C_azimuth_samples-NN+Y_c -1, abs(point_target_c_of_s_ac));
title('point target C image power');
xlabel('Range Time domain(Samples)');
ylabel('Azimuth Time domain(Samples)');  
% --------------------------------------------------------------------------------------------------------