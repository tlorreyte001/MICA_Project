function [ data2,R2,Q,S ] = Pan_and_Tompkins( data, Fs )
%% Apply Pan and Tompkins algorythm on data signal
    

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


%% Band-pass

% Filter
filtered_signal_bandpass = filter(num_high, den_high, filtered_signal_low); 
% TF_filtered_signal_bandpass = TF_filtered_signal_low .* TF_filtered_signal_high;


%% Derivative

% Filter
Ts = 1/Fs;
num_deriv = [1 2 0 -2 1];
den_deriv = 8*Ts;
 
tf(num_deriv,den_deriv);

% Signal
differentiated_signal = filter(num_deriv, den_deriv, filtered_signal_bandpass);
% TF_differentiated_signal = fftshift(abs(fft(differentiated_signal,nfft)));


%% Square

% Signal
signal_square = abs(differentiated_signal).^2;

%% MWI
% Convolution with a square signal  

% Signal 
M = floor(0.1*Fs);
porte = ones(1,M);
SMW = conv(porte,signal_square)/M;

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

data2=[data2 0 0 0];

end

