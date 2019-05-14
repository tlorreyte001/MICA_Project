function [ result ] = Fibrillation( R2 )
%FIBRILATION

delta_RR=[];
for i=1:length(R2)-1
    delta_RR = [delta_RR R2(i+1)-R2(i)];
end

delta_barre = mean(delta_RR);
gamma_estim = [];

for k=1:length(delta_RR)
    A=[];
    for n = 1:(length(delta_RR)-k)
        A=[A (delta_RR(n+k)-delta_barre)*(delta_RR(n)-delta_barre)];
    end
    gamma_estim(k) = (1/(length(delta_RR)-k))*sum(A);
end

if sum(gamma_estim>5000)>=1
    result = 'Atrial fibrillations';
else
    result = 'No fibrillation';
end

end

