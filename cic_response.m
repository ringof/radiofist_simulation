%% Copyright (c) 2015 by David Goncalves <davegoncalves@gmail.com>
%% See LICENCE.txt for details
%%
%% Script to calculate and plot the frequency. phase and impuslse response of a CIC
%% filter given the filter design parameters
%% See Altera Application Note AN455, "Understanding CIC Compensation Filters"
 
%% These are the design parameters for the CIC filter

M = 1;  	%% Differential Delay
N = 4; 		%% Number of Integrator/Comb Filter Pairs
R = 16; 	%% Sample Rate Decimation Ratio

w = 0:1e-3:pi/4;
z = e.^(1j * w);

% CIC filter frequency response

H = (1/R * (1 - z .^ (-R*M)) ./ (1 - z .^ (-1))) .^ N;

% Plot CIC filter frequency response

H_mag = 20 * log10(abs(H));
H_phase = atan2(imag(H),real(H));
%% TODO: Impulse response

%% Magnitude and Phase Plots

figure

subplot(2,1,1)
plot(w, H_mag)
title('CIC Filter Magnitude Response')
ylabel('|H| in db')
xlabel('w in Hz')
axis([min(w) max(w) -100 0])
grid on

subplot(2,1,2)
plot(w, H_phase)
title('CIC Filter Phase Response')
ylabel('angle in rad')
xlabel('w in Hz')
axis([min(w) max(w) -pi pi])
grid on

figure
plot(C_imp)
title('CIC Filter Impulse Response')
ylabel('amplitude')
xlabel('time')
grid on
