wavelengths = [800 840 880 920 960 1000 1040 1080 1120 1160 1200 1225 1250 1275 1300];
wavelengths = [1100 1110 1120 1130 1140 1150 1160 1170 1180]; 
wavelengths = [1120 1130]; 
wavelengths = [1120 1130 1140 1150 1160 1170 1180 1190 1200 1210 1215 1220 1225 1230 1240 1250 1260 1270 1280 1290 1300]; 
% wavelengths = [1210 1220 1230 1240 1250];
brightnesses = zeros(size(wavelengths)); 
bleachrates = brightnesses; 
laserpower = brightnesses; 
[background, ~] = jfbleach_analyze(0, 0, 4);

for i = 1:length(wavelengths)
	[br, ble] = jfbleach_analyze(wavelengths(i), background, 4); 
	brightnesses(i) = br; 
	bleachrates(i) = ble;  
end

% ran the experiment at 11% AOM, 
% measured power at 100%. 
power_meas = [45
				50
				72.5
				80.2
				69.2
				61.5
				52.2
				35.7
				24
				16.9
				8.8
				6.25
				4.7
				3.1
				1.5];
			
% run 4. 
% 69.1 @ 940 -- set the AOM at 50% to get 30.5mw at the sample.
% 1120 1130 1140 1150 1160 1180 1200 1220 1240 1260 1280 1300
% 100% AOM
% 1120 1130 1140 1150 1160 1170 1180 1190 1200 1210 1215 
% 1220 1225 1230 1240 1250 1260 1270 1280
power_meas = [
	42.8 % 1120
	43.8 % 1130
	44.8 % 1140
	44.4 % 1150
	44.5 % 1160
	45.8 % 1170
	47.1 % 1180
	41.3 % 1190
	35.4 % 1200
	33.2 % 1210
	32.1 % 1215
	31.0 % 1220
	30.9 % 1225
	30.8 % 1230
	30.9 % 1240
	30.7 % 1250
	30.5 % 1260
	26.9 % 1270
	23.3 % 1280
	19.6 % 1290
	15.9]; % 1300

% 940nm: 30% = 18mW at the sample
% 1230nm: 62% = 18mW at the sample
% 940nm: 47% = 30mW at the sample
% 1230nm: 100% = 30.4mW at the sample
% 1040nm, fixed output: 29% = 30.5mW at the sample. 
% 1040nm, variable: 60% = 30.4mW at the sample. 
			
subplot(1, 2, 1)
hold on
% plot(wavelengths, brightnesses ./ power_meas', 'Linewidth', 2)
plot(wavelengths, brightnesses, 'r', 'Linewidth', 2)
title('brightness, photon counts / mW excitation'); 
subplot(1,2,2)
hold on
plot(wavelengths, bleachrates, 'r', 'Linewidth', 2)
title('bleach rate, %/minute (approx)')