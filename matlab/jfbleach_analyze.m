function [brightness, bleach_rate, laser_power] = ...
	jfbleach_analyze(wavelength, background, ver)
% bleach rate is % / minute. 
load(['jf669_bleach' num2str(ver) '_' num2str(wavelength) '.mat']); 

brightness = mean(dat(:, 1)) - background; 

d2 = (dat(:,1)-background) / brightness; 

time = (dat(:, 2) - dat(1,2)) / 60; 

A = [time ones(size(dat, 1), 1)]; 
B = d2; 
coef = A\B; 
bleach_rate = -100 * coef(1); 
