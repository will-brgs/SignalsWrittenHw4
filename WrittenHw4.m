%% Signals and Systems Written Homework #4
%% Introduction
% * Author:                   Will Burgess
% * Class:                    ESE 351
% * Date:                     Created 2/17/2024, Last Edited 2/26/2024
% * With contributions from:  Mack Larosa, Tasha Igic, Mischa Tranor
%%
%% a
N = (2*pi) / (pi/5);
n = 0:1:N-1; 
a = 3 + sin(4*pi/5*n + pi/10) + cos(2*pi*n) + (-1).^n;

afft = fft(a);

ak_a = afft/N;

figure;
hold on
subplot(2,1,1);
stem(abs(fftshift(afft)), LineWidth=1.5);
title('Real Part')
xlabel('k value');
ylabel('Ak Output');
subplot(2,1,2)
stem(angle(fftshift(afft)), LineWidth=1.5);
title('Imaginary Part')
xlabel('k value');
ylabel('Ak Output');
sgtitle('FFT For 1.a)');

hold off
%% b, one period
b = [1,0,-2,0,0];
N = 5;
bfft = fft(b);
ak_b = bfft / N;

figure;
hold on
subplot(2,1,1);
stem(real(fftshift(b)), LineWidth=1.5);
title('Real Part')
xlabel('k value');
ylabel('Ak Output');
subplot(2,1,2)
stem(imag(fftshift(b)), LineWidth=1.5);
title('Imaginary Part')
xlabel('k value');
ylabel('Ak Output');
sgtitle('FFT For 1.b)');
hold off

%% C, N=4
N = 4;
n = 0:1:N-1;
c1 = 1 - sin(pi/2*n);
c1fft = fft(c1);
ak_c1 = c1fft/N;

figure;
hold on
subplot(2,1,1);
stem(real(fftshift(c1)), LineWidth=1.5);
title('Real Part')
xlabel('k value');
ylabel('Ak Output');
subplot(2,1,2)
stem(imag(fftshift(c1)), LineWidth=1.5);
title('Imaginary Part')
xlabel('k value');
ylabel('Ak Output');
sgtitle('FFT For 1.C: N = 4');
hold off

%% C, N=16
N = 16;
n = 0:1:N-1;
c2 = 1 - sin(pi/2*n);
c2fft = fft(c2);
ak_c2 = c2fft/N;

figure;
hold on
subplot(2,1,1);
stem(real(fftshift(c2)), LineWidth=1.5);
title('Real Part')
xlabel('k value');
ylabel('Ak Output');
subplot(2,1,2)
stem(imag(fftshift(c2)), LineWidth=1.5);
title('Imaginary Part')
xlabel('k value');
ylabel('Ak Output');
sgtitle('FFT For 1.C: N = 14');
hold off

%% d
N = (2*pi)/(pi/4);
n = 0:1:N-1; 
d = sin(7*pi/2*n) + exp(1j*pi/4*n);
dfft = fft(d);
ak_d = dfft/N;

figure;
hold on
subplot(2,1,1);
stem(real(fftshift(d)), LineWidth=1.5);
title('Real Part')
xlabel('k value');
ylabel('Ak Output');
subplot(2,1,2)
stem(imag(fftshift(d)), LineWidth=1.5);
title('Imaginary Part')
xlabel('k value');
ylabel('Ak Output');
sgtitle('FFT For 1.d');
hold off

%% Question 2:
N = 8;
%% a
k = -4:1:3;
ak = [-1,-1j,0,3,2,3,0,1j];

aifft = ifft(ak);
hand = 

figure
hold on
subplot (2,1,1)
stem(k,(fftshift(aifft), LineWidth=1.5);
title('IFFT for a');
xlabel('k value');
ylabel('Ak Output');
subplot (2,1,2)
stem(k,hand, LineWidth=1.5);
title('Hand Calcualation for b');
xlabel('k value');
ylabel('Ak Output');
sgtitle('IFFT for 2.a')
hold off
%% b
k = 0:1:N;
ak = cos((pi*k)/4);

bifft = ifft(ak);

figure
hold on
subplot (2,1,1)
stem(k,ifftshift(bifft), LineWidth=1.5);
title('IFFT for a');
xlabel('k value');
ylabel('Ak Output');
subplot (2,1,2)
stem(k,hand, LineWidth=1.5);
title('Hand Calcualation for c');
xlabel('k value');
ylabel('Function Output');
sgtitle('IFFT for 2.b')
hold off
%% c
k = -2:1:6;
ak = [1,1,1,1,1,0,0,0];
cbifft = ifft(ak);

subplot(3,1,3);
figure
hold on
subplot (2,1,1)
stem(k,real(ifftshift(cifft)), LineWidth=1.5);
title('IFFT for a');
xlabel('k value');
ylabel('Ak Output');
subplot (2,1,2)
stem(k,hand, LineWidth=1.5);
title('IFFT for a');
xlabel('k value');
ylabel('Function Output');
sgtitle('IFFT for 2.c')
hold off
