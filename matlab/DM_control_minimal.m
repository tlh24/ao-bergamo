% minimal closed-loop control of the DM

load('../data/calibration_flat.mat', 'mask');
load('../data/calibration_forward.mat', 'Cforward', 'cmask'); % valid lenslets in the actual microscope
load('../data/calibration_geneopt.mat', 'avg_wfs_x', 'avg_wfs_y');
load('../data/calibration_zernike.mat'); % Z dZx dZy 

nlenslets = sum(mask); 

calib = [avg_wfs_x avg_wfs_y];

DMcommand = zeros(97, 1); 

mmf = memmapfile('../shared_centroids.dat','Format','single','Offset',0,'Repeat',6000);

dmctrl = memmapfile('../shared_dmctrl.dat','Format','single','Offset',0,'Repeat',97, 'Writable',true);

% get the centroid positions. 
dat = mmf.Data; 
datx = dat(1:2:end); 
daty = dat(2:2:end); 

dx = datx(cmask) - calib(:,1); 
dy = daty(cmask) - calib(:,2);

des_coef = zeros(size(Z, 2), 1); 
des_coef(5) = 0.5*sin(n / 10); % spherical
des_dx = dZx * des_coef; 
des_dy = dZy * des_coef; 

A = [(dx - mean(dx)) - des_dx; ...
		 (dy - mean(dy)) - des_dy; 1]; %remove tip/tilt
cmd = A' * Cforward ; 
cmd = cmd' * -0.2; 
DMcommand = 0.999*DMcommand + cmd; 
% don't need to remove tip/tilt; shwfs already has this. 
% shwfs also clips to +- 0.15

% dmctrl.Data = single(DMcommandP); 