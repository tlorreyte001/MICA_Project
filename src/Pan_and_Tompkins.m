%% Mehdi GATI ; Thomas LORREYTE
clear; 
close all; clc;
addpath(genpath('.'));

%% Load a signal
[file,path] = uigetfile('*.mat', 'rt');
signal = load(fullfile(path, file));
data = signal.ecg; % Your ecg data
Fs = signal.Fs; % Sampling frequency
N = size(data,2); % Data length
time_axis = (1:N)/Fs; %% Signal
nfft = N; %fft points
axe_freq = linspace(-Fs/2,Fs/2, nfft);

data=-data; 

% Display
% figure;
% plot(time_axis, data); grid on;
% % hold on; plot(time_axis, th*ones(1,N), 'red');
% xlabel('Time (s)');
% ylabel('Magnitude');
% title('Time evolution of the loaded signal')

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
% subplot(3,1,1);
% plot(axe_freq, fftshift(abs(fft(data))));
% subplot(3,1,2);
% plot(axe_freq, TF_filtered_signal_low);
% subplot(3,1,3);
% plot(axe_freq, TF_filtered_signal_high);

%% Band-pass

% Filter
filtered_signal_bandpass = filter(num_high, den_high, filtered_signal_low); 
TF_filtered_signal_bandpass = TF_filtered_signal_low .* TF_filtered_signal_high;

% Display
% figure;
% subplot(2,1,1);
% plot(axe_freq, fftshift(abs(fft(data))));
% subplot(2,1,2);
% plot(axe_freq, TF_filtered_signal_bandpass);

%% Derivative

% Filter
Ts = 1/Fs;
num_deriv = [1 2 0 -2 1];
den_deriv = 8*Ts;

% Signal
differentiated_signal = filter(num_deriv, den_deriv, filtered_signal_bandpass);
TF_differentiated_signal = fftshift(abs(fft(differentiated_signal,nfft)));

% Display
% figure;
% subplot(2,1,1);
% plot(axe_freq, TF_filtered_signal_bandpass);
% subplot(2,1,2);
% plot(axe_freq, TF_differentiated_signal);

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
% figure;
% plot(SMW(23:end)/mean(abs(SMW(23:end))),'blue');
% hold on
% plot(data(1:length(SMW(23:end)))/max(abs(data)),'red');
% hold on


%% R-Searching 

SMW=SMW(23:end)/mean(abs(SMW)); % SMW without delay and "normalized" with the average value

%Searching intervals of rising slopes
[blue_peaks_amp, blue_peaks_loc] = findpeaks(SMW);
slope1=[];
slope1_amp=[];
slope2=[];
slope2_amp=[];
for i=1:length(blue_peaks_loc)-1
    if blue_peaks_amp(i+1)-blue_peaks_amp(i)>2 %Threshold for SMW rising slopes 
        slope1_amp=[slope1_amp blue_peaks_amp(i)]; 
        slope1=[slope1 blue_peaks_loc(i)];
        slope2=[slope2 blue_peaks_loc(i+1)]; 
        slope2_amp=[slope2_amp blue_peaks_amp(i+1)]; 
    end 
end

% Display
% plot(slope1,slope1_amp,'o');
% hold on;
% plot(slope2,slope2_amp,'o');

%R-peaks
data2 = data(1:length(SMW))/max(abs(data)); %Data normalized

data_slope_interval=[];
for i = 1:length(slope1)
    data_slope_interval=[data_slope_interval data2(slope1(i):slope2(i))];
end
red_treshold = mean(abs(data_slope_interval)); %Threshold for data rising peaks

R2 = [];
R2_amp=[];
for i=1:length(slope1)
    data_slope_interval=[data2(slope1(i):slope2(i))];
    [red_peaks_amp, red_peaks_loc] = findpeaks(data2(slope1(i):slope2(i)),'MinPeakHeight',red_treshold);
    if (~isempty(red_peaks_loc))
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
elseif (1/delta_barre*Fs*60)>100
    result='Tachycardia';
else
    result='Sain';
end
result

% Ectopic beat
ectopic = [];
for i=1:length(delta_RR)-1
    if delta_RR(i) < delta_barre-2*sqrt(var(delta_RR))
        ectopic=[ectopic R2(i)];
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
    gamma_estim(k) = (1/(length(delta_RR)-k))*sum(A);
end
gamma_estim;
figure;
plot(gamma_estim)