function [ result ] = Fibrilation( delta_RR )
%FIBRILATION Summary of this function goes here
%   Detailed explanation goes here

gamma_estim = [];
delta_barre = mean(delta_RR);

for k=1:length(delta_RR)
    A=[];
    for n = 1:(length(delta_RR)-k)
        A=[A (delta_RR(n+k)-delta_barre)*(delta_RR(n)-delta_barre)];
    end
    gamma_estim(k) = (1/(length(delta_RR)-k))*sum(A);
end

if sum(gamma_estim>5000)>=1
    result = 'Atrial Fibrilations';
else
    result = 'No fibrilation';
end

end

