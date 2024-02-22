%% Signals and Systems Written Homework #4
%% Introduction
% * Author:                   Will Burgess
% * Class:                    ESE 351
% * Date:                     Created 2/17/2024, Last Edited 2/26/2024
% * With contributions from:  Mack Larosa, Tasha Igic, Mischa Tranor
%%

% Define signal a
fs=100;
t = 0:1/fs:2*pi; 
a = 3 + sin(4*pi/5*t + pi/10) + cos(2*pi*t) + (-1).^t;

% Define signal b (delta function)
N = 50; % Number of periods
b = zeros(1, N);
b(1:5:end) = 1;
b(2:5:end) = -2;

% Define signal c (i) with period N=4
N_c1 = 4;
t_c1 = 0:1/fs:N_c1-1/fs;
c1 = 1 - sin(pi/2*t_c1);

% Define signal c (ii) with period N=16
N_c2 = 16;
t_c2 = 0:1/fs:N_c2-1/fs;
c2 = 1 - sin(pi/2*t_c2);

% Define signal d
d = sin(7*pi/2*t) + exp(1j*pi/4*t);

% Calculate fundamental frequency and coefficients for signal a
ca = fftshift(fft(a, N));
[~, I] = max(abs(ca));
fa = fs * (0:N-1) / N;  % Frequency vector

% Calculate fundamental frequency and coefficients for signal b
[fb, cb] = fftshift(fft(b, N));
[~, I] = max(abs(cb));
fb = fs * (I-1) / N;

% Calculate fundamental frequency and coefficients for signal c (i)
[fc1, cc1] = fftshift(fft(c1, N_c1));
[~, I] = max(abs(cc1));
fc1 = fs * (I-1) / N_c1;

% Calculate fundamental frequency and coefficients for signal c (ii)
[fc2, cc2] = fftshift(fft(c2, N_c2));
[~, I] = max(abs(cc2));
fc2 = fs * (I-1) / N_c2;

% Calculate fundamental frequency and coefficients for signal d
[fd, cd] = fftshift(fft(d, N));
[~, I] = max(abs(cd));
fd = fs * (I-1) / N;

% Plot signal a and its magnitude spectrum
figure(1);
subplot(2, 1, 1);
plot(t, a);
title('Signal a');
subplot(2, 1, 2);
plot(fa/fs, abs(ca));
title('Magnitude Spectrum of Signal a');

% Plot signal b and its magnitude spectrum
figure(2);
subplot(2, 1, 1);
plot(t(1:N), b);
title('Signal b');
subplot(2, 1, 2);
plot(fb/fs, abs(cb));
title('Magnitude Spectrum of Signal b');

% Plot signal c (i) and its magnitude spectrum
figure(3);
subplot(2, 1, 1);
plot(t_c1, c1);
title('Signal c (N=4)');
subplot(2, 1, 2);
plot(fc1/fs, abs(cc1));
title('Magnitude Spectrum of Signal c (N=4)');

% Plot signal c (ii) and its magnitude spectrum
figure(4);
subplot(2, 1, 1);
plot(t_c2, c2);
title('Signal c (N=16)');
subplot(2, 1, 2);
plot(fc2/fs, abs(cc2));
title('Magnitude Spectrum of Signal c (N=16)');

% Plot signal d and its magnitude spectrum
figure(5);
subplot(2, 1, 1);
plot(t, d);
title('Signal d');
subplot(2, 1, 2);
plot(fd/fs, abs(cd));

