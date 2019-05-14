function [ filtered_signal_bandpass, TF_filtered_signal_bandpass ] = band_pass( data )
%% BAND_PASS FILTER

%% High and low filters

% Filters

% Low pass filter
num_low = [1 0 0 0 0 0 -2 0 0 0 0 0 1];
den_low = [1 -2 1];

% High pass filter
num_high = zeros(1,33);
num_high(1) = -1;
num_high(17) = 32;
num_high(18) = -32;
num_high(33) = +1;
den_high = [1 -1];

% Signals
filtered_signal_low = filter(num_low, den_low, data);
filtered_signal_high = filter(num_high, den_high, data);
TF_filtered_signal_low = fftshift(abs(fft(filtered_signal_low)));
TF_filtered_signal_high = fftshift(abs(fft(filtered_signal_high)));

%% Band-pass

% Filter
filtered_signal_bandpass = filter(num_high, den_high, filtered_signal_low); 
TF_filtered_signal_bandpass = TF_filtered_signal_low .* TF_filtered_signal_high;


end

