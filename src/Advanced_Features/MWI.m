function [ SMW ] = MWI( differentiated_signal )
%% MWI filter

% Square
% Signal
signal_square = abs(differentiated_signal).^2;

% MWI
% Signal 
M = floor(0.1*Fs);
porte = ones(1,M);
SMW = conv(porte,signal_square)/M;

end

