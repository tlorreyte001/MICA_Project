function [ result, bpm ] = Brady_Tachy( R2, Fs )
%% Detect Brady or Tachycardia with R-pics position

delta_RR=[];
for i=1:length(R2)-1
    delta_RR = [delta_RR R2(i+1)-R2(i)];
end
delta_barre = mean(delta_RR);

bpm = 1/delta_barre*Fs*60;

if (1/delta_barre*Fs*60)<60
    result='Bradycardia';
elseif (1/delta_barre*Fs*60)>100
    result='Tachycardia';
else
    result='Healthy';
end


end

