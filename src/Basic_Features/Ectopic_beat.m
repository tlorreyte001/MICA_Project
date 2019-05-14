function [ ectopic ] = Ectopic_beat( R2, Fs )
% Detect position of Ectopic beat.

delta_RR=[];
for i=1:length(R2)-1
    delta_RR = [delta_RR R2(i+1)-R2(i)];
end
delta_barre = mean(delta_RR);

ectopic = [];
for i=1:length(delta_RR)-1
    if delta_RR(i) < delta_barre-2*sqrt(var(delta_RR))
        ectopic=[ectopic R2(i)];
    end
end


end