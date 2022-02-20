% this file regenerates 
% ../data/calibration_zernike.mat
% based on geneopt.

load('../data/calibration_flat.mat');
load('../data/calibration_forward.mat'); % valid lenslets in the actual microscope
load('../data/calibration_960nm_20220217_geneopt.mat');
nlenslets = sum(mask);

cx = mean(genecalib(:, 1)); 
cy = mean(genecalib(:, 2)); 
calibx = genecalib(:, 1); 
caliby = genecalib(:, 2); 

calibx_n = calibx - cx; 
caliby_n = caliby - cy; 
calibx_n = (calibx_n / 727.0); 
caliby_n = (caliby_n / 727.0); 
nlenslets_n = sum(mask); 
[Z,dZx,dZy] = zernikes_and_derivatives_cartesian_OSA (calibx_n, caliby_n, 7, "zeros");
Z = squeeze(Z); 
dZx = squeeze(dZx); 
dZy = squeeze(dZy); 
save('../data/calibration_zernike.mat', 'Z', 'dZx', 'dZy'); 

% check
figure; 
plot(calibx_n.^2 + caliby_n.^2)
title('values should not exceed 1')
% this one got me in the past! 