%% Main script to test ecg function without gui
% This file computes a simple analysis of an ecg signal. You can use it to test the different processing methods. 
% This first version will plot the temporal signal, compute its cardiac rythma and display the different P, Q, R, S, T points for a specific segment.  

clear; close all; clc;
addpath(genpath('.'));

%% Load a signal
[file,path] = uigetfile('*.mat', 'rt');
signal = load(fullfile(path, file));
data = signal.ecg; % Your ecg data
Fs = signal.Fs; % Sampling frequency
N = size(data,2); % Data length
time_axis = (1:N)/Fs; 

%% Threshold method
% th = 200; % threshold
% i_seg = 10; % Segment number to plot
% 
% % Time plot
figure;
plot(time_axis, data); grid on;
% hold on; plot(time_axis, th*ones(1,N), 'red');
xlabel('Time (s)');
ylabel('Magnitude');
title('Time evolution of the loaded signal')
% 
% % Print BPM
% [bpm, R_locs] = bpm_threshold(data, th, Fs);
% % Figures PQRST
% [segment, P_loc, Q_loc, R_loc, S_loc, T_loc] = ecg_threshold(data, R_locs, i_seg);
% time_segment = (1:length(segment))/Fs;
% 
% figure;
% h = plot(time_segment, segment); grid on;
% hold on;
% plot(time_segment(P_loc),segment(P_loc), '*','Color','red'); text(time_segment(P_loc),segment(P_loc),' P ','Color','red','FontSize',14);
% plot(time_segment(Q_loc),segment(Q_loc), '*','Color','red'); text(time_segment(Q_loc),segment(Q_loc),' Q ','Color','red','FontSize',14);
% plot(time_segment(R_loc),segment(R_loc), '*','Color','red'); text(time_segment(R_loc),segment(R_loc),' R ','Color','red','FontSize',14);
% plot(time_segment(S_loc),segment(S_loc), '*','Color','red'); text(time_segment(S_loc),segment(S_loc),' S ','Color','red','FontSize',14);
% plot(time_segment(T_loc),segment(T_loc), '*','Color','red'); text(time_segment(T_loc),segment(T_loc),' T ','Color','red','FontSize',14);
% hold off;
% xlabel('Time (s)');
% ylabel('Magnitude');
% title('ECG segment characteristic')

%% Your turn : My new method !

%% Negative signal

i=1;
while i<length(data)
    if data(i)<-0.4
        data=-data;
        break
    end
    i=i+1;
end
    

%% High and low filters

% Filters

% Low pass filter
num_low = [1 0 0 0 0 0 -2 0 0 0 0 0 1];
den_low = [1 -2 1];
%fvtool(num_low, den_low);
% High pass filter
num_high = zeros(1,33);
num_high(1) = -1;
num_high(17) = 32;
num_high(18) = -32;
num_high(33) = +1;
den_high = [1 -1];
%fvtool(num_high, den_high);

% Signals
filtered_signal_low = filter(num_low, den_low, data);
filtered_signal_high = filter(num_high, den_high, data);
TF_filtered_signal_low = fftshift(abs(fft(filtered_signal_low)));
TF_filtered_signal_high = fftshift(abs(fft(filtered_signal_high)));

% Display

% plot(time_axis, filtered_signal)

nfft = N;
axe_freq = linspace(-1000,1000, nfft);
subplot(3,1,1);
plot(axe_freq, fftshift(abs(fft(data))));
subplot(3,1,2);
plot(axe_freq, TF_filtered_signal_low);
subplot(3,1,3);
plot(axe_freq, TF_filtered_signal_high);

%% Band-pass

% Filter
filtered_signal_bandpass = filter(num_high, den_high, filtered_signal_low); 
TF_filtered_signal_bandpass = TF_filtered_signal_low .* TF_filtered_signal_high;

% Display
figure;
subplot(2,1,1);
plot(axe_freq, fftshift(abs(fft(data))));
subplot(2,1,2);
plot(axe_freq, TF_filtered_signal_bandpass);

%% Derivative

% Filter
Ts = 1/Fs;
num_deriv = [1 2 0 -2 1];
den_deriv = 8*Ts;
 
tf(num_deriv,den_deriv)
% Signal
differentiated_signal = filter(num_deriv, den_deriv, filtered_signal_bandpass);
TF_differentiated_signal = fftshift(abs(fft(differentiated_signal,nfft)));

% Display
figure;
subplot(2,1,1);
plot(axe_freq, TF_filtered_signal_bandpass);
subplot(2,1,2);
plot(axe_freq, TF_differentiated_signal);

%% Square

% Signal
signal_square = abs(differentiated_signal).^2;

%% MWI
% Convolution with a square signal  

% Signal 
M = floor(0.1*Fs);
porte = ones(1,M);
SMW = conv(porte,signal_square)/M;

% Display
figure;
plot(SMW(23:end)/max(abs(SMW(23:end))),'blue');
hold on
plot(data(1:length(SMW(23:end)))/max(abs(data)),'red');
hold on


%% R-Searching 

[blue_peaks_amp, blue_peaks_loc] = findpeaks(SMW(23:end)/max(abs(SMW(23:end))));
slope1=[];
slope1_amp=[];
slope2=[];
slope2_amp=[];
for i=1:length(blue_peaks_loc)-1
    if blue_peaks_amp(i+1)-blue_peaks_amp(i)>0.2
        slope1_amp=[slope1_amp blue_peaks_amp(i)]; 
        slope1=[slope1 blue_peaks_loc(i)];
        slope2=[slope2 blue_peaks_loc(i+1)]; 
        slope2_amp=[slope2_amp blue_peaks_amp(i+1)]; 
    end 
end

plot(slope1,slope1_amp,'o');
hold on;
plot(slope2,slope2_amp,'o');

data2 = data(1:length(SMW(23:end)))/max(abs(data));

R2 = [];
R2_amp=[];
for i=1:length(slope1)
    [red_peaks_amp, red_peaks_loc] = findpeaks(data2(slope1(i):slope2(i)),'MinPeakHeight',0.15);
    if (length(red_peaks_loc)==0)
        R2_amp = [R2_amp 0];
        R2 = [R2 slope2(i)-1];
    else
        R2_amp = [R2_amp red_peaks_amp];
        R2=[R2 red_peaks_loc+slope1(i)-1];
    end
end

plot(R2,R2_amp,'o')
sum(R2_amp==0);
        
%% Q and S detection

Q=[];
for i=1:length(R2)
    j=1;
    while(data2(R2(i)-j)<data2(R2(i)-j+1))
        j=j+1;
    end
    Q=[Q R2(i)-j+1];
end


S=[];
for i=1:length(R2)
    j=1;
    while((data2(R2(i)+j)-data2(R2(i)+j+1))>0.01)
        j=j+1;
    end
    S=[S R2(i)+j];
end


figure;
plot(data2);
hold on;
plot(R2,R2_amp,'o');
plot(Q,data2(Q),'o');
plot(S,data2(S),'o');

%% Identification of cardiac pathologies

% Bradycardia/Tachycardia
delta_RR=[];
for i=1:length(R2)-1
    delta_RR = [delta_RR R2(i+1)-R2(i)];
end
delta_barre = mean(delta_RR);

if (1/delta_barre*Fs*60)<60
    result='Bradycardia';
elseif (delta_barre/Fs*60)>100
    result='Tachycardia';
else
    result='Sain';
end
result

% Ectopic beat
ectopic = [];
for i=1:length(delta_RR)-1
    if delta_RR(i) < delta_barre-2*sqrt(var(delta_RR))
        ectopic=[ectopic i];
    end
end
ectopic;

% Fibrillation

gamma_estim = [];

for k=1:length(delta_RR)
    A=[];
    for n = 1:(length(delta_RR)-k)
        A=[A (delta_RR(n+k)-delta_barre)*(delta_RR(n)-delta_barre)];
    end
    gamma_estim=[gamma_estim (1/(length(delta_RR)-k))*sum(A)];
end
gamma_estim;
        