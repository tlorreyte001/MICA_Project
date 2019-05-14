function [ SMW, data2 ] = MWI( differentiated_signal, data, Fs )
%% MWI filter

% Square
% Signal
signal_square = abs(differentiated_signal).^2;

% MWI
% Signal 
M = floor(0.1*Fs);
porte = ones(1,M);
SMW = conv(porte,signal_square)/M;

SMW=SMW(23:end)/max(abs(SMW));
data2 = data(1:length(SMW(23:end)))/max(abs(data));
data2 = [data2 0 0 0];

end

