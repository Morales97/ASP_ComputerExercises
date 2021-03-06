addpath(pwd)
cd 'Support scripts'/
addpath(pwd)
cd ..
cd 'Wiener filters'/
addpath(pwd)
cd ..
cd mfiles/
addpath(pwd)
cd ..
[z,fs] = audioread('EQ2401project1data2021.wav');


%% AR model order selection - signal

[z,fs] = audioread('EQ2401project1data2021.wav');
M_lim = 200;

sig = z;
mdl_aic(sig, M_lim);

%% AR model order selection - noise

[z,fs] = audioread('EQ2401project1data2021.wav');
x = z(1:8000);      % Noise samples
M_lim = 200;

sig = x;
mdl_aic(sig, M_lim);

%% FIR filter
x = z(1:8000);      % Noise samples
N = 30;             % FIR filter length
M_signal = 30;      % AR model order for z(n)
M_noise = 10;       % AR model order for x(n)

[shatfir, thetahatfir] = p1_firw(z, x, N, M_signal, M_noise);

%% Causal filter
x = z(1:8000);      % Noise samples
M_signal = 30;      % AR model order for z(n)
M_noise = 10;       % AR model order for x(n)

[shatc, numc, denc] = p1_cw(z, x, M_signal, M_noise);

%% Non-Causal filter
x = z(1:8000);      % Noise samples
N = 40;             % FIR filter length
M_signal = 30;      % AR model order for z(n)
M_noise = 10;       % AR model order for x(n)
BT_lag = 60;        % Blackman-Tuckey lag
use_BT = 0;         % '1' to use Blackman-tuckey PSD estimation, '0' to use 
                    % AR estimation

[shatnc, numnc, dennc] = p1_ncw(z, x, M_signal, M_noise, BT_lag, use_BT);
