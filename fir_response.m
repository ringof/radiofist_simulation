%% Copyright (c) 2015 by David Goncalves <davegoncalves@gmail.com>
%% See LICENCE.txt for details
%%
%% Script to calculate and plot the frequency. phase and impuslse response of a FIR
%% filter given the filter design parameters
%% This code is primarily taken from the Altera Application Note AP455, 
%% "Understanding CIC Compensation Filters", with some prettifying for my own
%% sanity. 
 
%% CIC Filter Design Parameters

M = 1;  	%% Differential Delay
N = 4; 		%% Number of Integrator/Comb Filter Pairs
R = 16; 	%% Sample Rate Decimation Ratio

%% FIR Filter Design Parameters

B = 18; 	%% Coeffi. Bit-width
Fs = 6e6; 	%% (High) Sampling freq in Hz before decimation
Fc = 5e3; 	%% Pass band edge in Hz
L = 100; %% Filter order; must be even
Fo = R*Fc/Fs; %% Normalized Cutoff freq; 0<Fo<=0.5/M;

%% CIC Compensator Design using fir2 
p = 2e3; %% Granularity
s = 0.25/p; %% Step size
fp = [0:s:Fo]; %% Pass band frequency samples
fs = (Fo+s):s:0.5; %% Stop band frequency samples
f = [fp fs]*2; %% Normalized frequency samples; 0<=f<=1

Mp = ones(1,length(fp)); %% Pass band response; Mp(1)=1
Mp(2:end) = abs( M*R*sin(pi*fp(2:end)/R)./sin(pi*M*fp(2:end))).^N;
Mf = [Mp zeros(1,length(fs))];
f(end) = 1;

h = fir2(L,f,Mf); %% Filter length L+1

[z, w] = freqz(h);
 
h = h/max(h); %% Floating point coefficients

hz = round(h*power(2,B-1)-1); %% Fixed point coefficients

% Plot FIR filter frequency response

H_mag = 20 * log10(abs(z));
H_phase = atan2(imag(z),real(z));
%% TODO: Impulse response

%% Magnitude and Phase Plots

figure

subplot(2,1,1)
plot(w/pi, H_mag)
title('Compensating FIR Filter Magnitude Response')
ylabel('|H| in db')
xlabel('w in Hz')
axis([min(w/pi) max(w/pi) -100 0])
grid on

subplot(2,1,2)
plot(w/pi, H_phase)
title('FIR Filter Phase Response')
ylabel('angle in rad')
xlabel('f/fs')
axis([min(w/pi) max(w/pi) -pi pi])
grid on

figure
plot(C_imp)
title('CIC Filter Impulse Response')
ylabel('amplitude')
xlabel('time')
grid on
