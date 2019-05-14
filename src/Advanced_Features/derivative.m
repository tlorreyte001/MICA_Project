function [ differentiated_signal, TF_differentiated_signal ] = derivative( filtered_signal_bandpass )
%% Derivative filter

% Filter
Ts = 1/Fs;
num_deriv = [1 2 0 -2 1];
den_deriv = 8*Ts;

% Signal
differentiated_signal = filter(num_deriv, den_deriv, filtered_signal_bandpass);
TF_differentiated_signal = fftshift(abs(fft(differentiated_signal,nfft)));
end

