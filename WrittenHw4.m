%% Signals and Systems Written Homework #4
%% Introduction
% * Author:                   Will Burgess
% * Class:                    ESE 351
% * Date:                     Created 2/17/2024, Last Edited 2/26/2024
% * With contributions from:  Mack Larosa, Tasha Igic, Mischa Tranor
close all
%% Question 1
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
stem(real(fftshift(bfft)), LineWidth=1.5);
title('Real Part')
xlabel('k value');
ylabel('Ak Output');
subplot(2,1,2)
stem(imag(fftshift(bfft)), LineWidth=1.5);
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
stem(real(fftshift(c1fft)), LineWidth=1.5);
title('Real Part')
xlabel('k value');
ylabel('Ak Output');
subplot(2,1,2)
stem(imag(fftshift(c1fft)), LineWidth=1.5);
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
stem(real(fftshift(c2fft)), LineWidth=1.5);
title('Real Part')
xlabel('k value');
ylabel('Ak Output');
subplot(2,1,2)
stem(imag(fftshift(c2fft)), LineWidth=1.5);
title('Imaginary Part')
xlabel('k value');
ylabel('Ak Output');
sgtitle('FFT For 1.C: N = 16');
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
stem(real(fftshift(dfft)), LineWidth=1.5);
title('Real Part')
xlabel('k value');
ylabel('Ak Output');
subplot(2,1,2)
stem(imag(fftshift(dfft)), LineWidth=1.5);
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

figure, hold on
subplot(3,2,1)
stem(k,real(ifftshift(aifft)),'b',LineWidth=1.5);
title('Real Part IFFT for a');
xlabel('k value');
ylabel('Function Output');

subplot(3,2,2)
stem(k,imag(ifftshift(aifft)), 'r',LineWidth=1.5);
title('Imaginary Part IFFT for a');
xlabel('k value');
ylabel('Function Output');

%% b
k = 0:1:N;
ak = cos((pi*k)/4);

bifft = ifft(ak);

subplot(3,2,3)
stem(k,real(ifftshift(bifft)), 'b',LineWidth=1.5);
title('Real Part IFFT for a');
xlabel('k value');
ylabel('Function Output');

subplot(3,2,4)
stem(k,imag(ifftshift(bifft)),'r', LineWidth=1.5);
title('Imaginary Part IFFT for a');
xlabel('k value');
ylabel('Function Output');

%% c
k = -2:1:6;
ak = [1,1,1,1,1,0,0,0,0];
cifft = ifft(ak);

subplot(3,2,5)
stem(k,real(ifftshift(cifft)), 'b',LineWidth=1.5);
title('Real Part IFFT for a');
xlabel('k value');
ylabel('Function Output');

subplot(3,2,6)
stem(k,imag(ifftshift(cifft)),'r', LineWidth=1.5);
title('Imaginary Part IFFT for a');
xlabel('k value');
ylabel('Function Output');
sgtitle('Queston 2 Function Outputs')
hold off
%% Question 3b
%% 3bi
N1 = 3;
N = [16,32,64];

% N = 16
n = (-N1:1:N(1)-N1);
x = zeros(1,length(n));
x(1:N(1)-N1)= 1;

ft = fft(x);

figure, hold on
subplot (3,2,1)
stem(n,abs(ft),'b', LineWidth=1.5);
title('Magnitude for N= 16');
xlabel('index n');
ylabel('Fourier Series Coefficient');
subplot (3,2,2)
stem(n,angle(ft), 'r', LineWidth=1.5);
title('Phase for N= 16');
xlabel('index n');
ylabel('Fourier Series Coefficient');

%N = 32
n = (-N1:1:N(2)-N1);
x = zeros(1,length(n));
x(1:N(2)-N1)= 1;

ft = fft(x);

subplot (3,2,3)
stem(n,abs(ft),'b', LineWidth=1.5);
title('Magnitude for N= 32');
xlabel('index n');
ylabel('Fourier Series Coefficient');
subplot (3,2,4)
stem(n,angle(ft), 'r', LineWidth=1.5);
title('Phase for N= 32');
xlabel('index n');
ylabel('Fourier Series Coefficient');

% N= 64
n = (-N1:1:N(3)-N1);
x = zeros(1,length(n));
x(1:N(3)-N1)= 1;

ft = fft(x);

subplot (3,2,5)
stem(n,abs(ft),'b', LineWidth=1.5);
title('Magnitude for N= 64');
xlabel('index n');
ylabel('Fourier Series Coefficient');
subplot (3,2,6)
stem(n,angle(ft), 'r', LineWidth=1.5);
title('Phase for N= 64');
xlabel('index n');
ylabel('Fourier Series Coefficient');
sgtitle('Fourier Series Coefficients for bi:')
hold off

% Observation : We see identical relationships with the fourier series
% coefficients with respective chagnes to N, notably, magnitude has a very
% large magnitude at the first value calculated, while the rest are low and
% have a roughly wave shaped appearence. For the phase we see an identical
% triangular wedge pattern which is given more datapoints in higher N1
% values


%% 3bii
N1vector = [2,6,10];
N = 32;

% N1 = 2
N1 = N1vector(1);
n = (-N1:1:N-N1);
x = zeros(1,length(n));
x(1:N-N1)= 1;

ft = fft(x);

figure, hold on
subplot (3,2,1)
stem(n,abs(ft),'b', LineWidth=1.5);
title('Magnitude for N= 2');
xlabel('index n');
ylabel('Fourier Series Coefficient');
subplot (3,2,2)
stem(n,angle(ft), 'r', LineWidth=1.5);
title('Phase for N= 2');
xlabel('index n');
ylabel('Fourier Series Coefficient');

%N1 = 6

N1 = N1vector(2);
n = (-N1:1:N-N1);
x = zeros(1,length(n));
x(1:N-N1)= 1;

ft = fft(x);

subplot (3,2,3)
stem(n,abs(ft),'b', LineWidth=1.5);
title('Magnitude for N1= 6');
xlabel('index n');
ylabel('Fourier Series Coefficient');
subplot (3,2,4)
stem(n,angle(ft), 'r', LineWidth=1.5);
title('Phase for N1= 6');
xlabel('index n');
ylabel('Fourier Series Coefficient');

%N1 = 10

N1 = N1vector(3);
n = (-N1:1:N-N1);
x = zeros(1,length(n));
x(1:N-N1)= 1;

ft = fft(x);

subplot (3,2,5)
stem(n,abs(ft),'b', LineWidth=1.5);
title('Magnitude for N1= 10');
xlabel('index n');
ylabel('Fourier Series Coefficient');
subplot (3,2,6)
stem(n,angle(ft), 'r', LineWidth=1.5);
title('Phase for N1= 10');
xlabel('index n');
ylabel('Fourier Series Coefficient');
sgtitle('Fourier Series Coefficients for bii:')

% Observation: While the output bears similarity to bi, there are clear
% differences. Magnitude shares the intiial spike but shows a concave
% upwards trend close to the enpoint of indicies run through fft. For phase
% we see the wedge shape seen before but as N1 increases there is large
% distortion. For N1 we cannot tell the wedge shape at all

%% 3biii
N1 = 3;
N = 16;
n_reg = (-N1:1:N-N1);
x_reg = zeros(1,length(n));
x_reg(1:N-N1)= 1;

% x[n-3]
n = (-N1-3:1:N-N1-3);
n = n - 3;
x = zeros(1,length(n));
x(1:N-N1) = 1;
ft = fft(x);

figure, hold on
subplot (3,2,1)
stem(n,abs(ft),'b', LineWidth=1.5);
title('Magnitude for x[n-3]');
xlabel('index n');
ylabel('Fourier Series Coefficient');
subplot (3,2,2)
stem(n,angle(ft), 'r', LineWidth=1.5);
title('Phase for x[n-3]');
xlabel('index n');
ylabel('Fourier Series Coefficient');

% x[n]-1/2
n = (-N1:1:N-N1);
x = zeros(1,length(n));
x(1:N-N1) = 1;
x = x - 1/2;
ft = fft(x);

subplot (3,2,3)
stem(n,abs(ft),'b', LineWidth=1.5);
title('Magnitude for x[n]-1/2');
xlabel('index n');
ylabel('Fourier Series Coefficient');
subplot (3,2,4)
stem(n,angle(ft), 'r', LineWidth=1.5);
title('Phase for x[n]-1/2');
xlabel('index n');
ylabel('Fourier Series Coefficient');

%cos(pi*n)*x[n]
n = (-N1:1:N-N1);
x = zeros(1,length(n));
x(1:N-N1) = 1;
x = x .* cos(pi*n);
ft = fft(x);

subplot (3,2,5)
stem(n,abs(ft),'b', LineWidth=1.5);
title('Magnitude for x[n]cos(pi*n)');
xlabel('index n');
ylabel('Fourier Series Coefficient');
subplot (3,2,6)
stem(n,angle(ft), 'r', LineWidth=1.5);
title('Phase for x[n]cos(pi*n)');
xlabel('index n');
ylabel('Fourier Series Coefficient');
sgtitle('Fourier Series Coefficients for biii')
hold off

%Observation: The wedge shapes described in earlier parts of 3b are also
%seen here. We can see that the shift for the first plot makes no change to
%the resulting fourier series coefficients. In fact, both the first and
%second plots are identical to each other, showing the transformations of
%x[n] have no bearing on the resulting ak. For the cosine function, I
%expected no change but was wrong. The wedge shape is harder to notice and
%the coefficients are seen to be symmetric as expected
%% Question 3c
%% 3ci

