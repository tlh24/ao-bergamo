%%--- 
load('../calibration_flat.mat');
load('../calibration_forward.mat'); % valid lenslets in the actual microscope

nlenslets = size(calib, 1); 
mask = mask(1:nlenslets); 
% select 'valid' centroids based on calibration data. 
cx = mean(calib(mask, 1)); 
cy = mean(calib(mask, 2)); 
calibx = calib(mask, 1); 
caliby = calib(mask, 2); 

calibx_n = calibx - cx; 
caliby_n = caliby - cy; 
calibx_n = (calibx_n / 925.0); 
caliby_n = (caliby_n / 925.0); 
nlenslets_n = sum(mask); 
[Z,dZx,dZy] = zernikes_and_derivatives_cartesian_OSA (calibx_n, caliby_n, 6, "NaN");

mmf = memmapfile('../shared_centroids.dat','Format','single','Offset',0,'Repeat',6000);

dmctrl = memmapfile('../shared_dmctrl.dat','Format','single','Offset',0,'Repeat',97, 'Writable',true);

% doctor the forward control matrix to remove centroids where the mean
% weight is more than a standard deviation off. 
Cmask = [goodmask; goodmask; 1]; 
Cmask = repmat(Cmask, 1, 97); 
Cforward = C .* Cmask; 
% also clamp the values... 
Cforward = min(Cforward, 0.03); 
Cforward = max(Cforward, -0.03); 

DMcommand = zeros(97, 1); 
[dmx, dmy] = dm_actuator_to_xy();

fig1 = figure;

for x = 1:300
    dat = mmf.Data; 
    dx = dat(1:2:end); 
    dy = dat(2:2:end); 
    
    dx = dx(1:nlenslets) - calib(:,1); 
    dy = dy(1:nlenslets) - calib(:,2); 
    
    dx = dx(mask); 
    dy = dy(mask); 

    % use these derivatives to calculate Zernike coefficients. 

    % serially regress dx and dy to the cartesian derivatives. 
    coef = zeros(28, 1); 
    B = [dx ; dy]; 
    indx = 2; 
    len = 1; 
    while(indx < 28)
        fin = indx + len; 
        % disp(['working on ' num2str(indx) ' to ' num2str(fin)]); 
        A = []; 
        for i = indx : fin
            A = [A [dZx(:,1,i); dZy(:,1,i)]];
        end
        C = A\B; 
        pred = A*C; 
        % visually check!
        if 0
            figure; 
            subplot(1,2,1); 
            scatter3(calibx_n, caliby_n, B(1:nlenslets_n)); 
            hold on
            scatter3(calibx_n, caliby_n, pred(1:nlenslets_n), 'ro'); 
            title(['dx fit for Zernike modes ' num2str(indx) ' to ' num2str(fin)]);
            subplot(1,2,2); 
            scatter3(calibx_n, caliby_n, B(nlenslets_n+1:end)); 
            hold on
            scatter3(calibx_n, caliby_n, pred(nlenslets_n+1:end), 'ro'); 
            title(['dy fit for Zernike modes ' num2str(indx) ' to ' num2str(fin)]);
        end
        B = B - pred; 
        coef(indx:fin) = C; 
        indx = fin+1;
        len = len+1;
    end

    % looks alright, I guess. 
    % try reconstructing based on the Zernike polynomials. 
    Z = squeeze(Z); 
	 coef2 = coef; 
	 coef2(2:3) = 0; 
    D = Z * coef2; 
    figure(fig1); 
    subplot(1,3,1); 
    scatter3(calibx_n, caliby_n, D/max(D)); 
    title('Reconstruction')
    subplot(1,3,2); 
    bar(coef); 
    title('Zernike coefficients'); 
    
    % now, based on dx and dy, calculate updated control signals .. 
    A = [dx - mean(dx); dy - mean(dy); 1]; %remove tip/tilt
    cmd = A' * Cforward ; 
	 cmd = cmd' * -0.025; 
 	 DMcommand = 0.95*DMcommand + cmd; 
	 DMcommand = min(DMcommand, 0.15); % clip, per reality.
	 DMcommand = max(DMcommand, -0.15); 
	 DMcommand = DMcommand - mean(DMcommand); %piston
	 tilt = dmx \ DMcommand;
	 % there still could be tip-tilt in the command, 
	 % even if it was removed from the WFS derivatives.
	 DMcommand = DMcommand - dmx * tilt; 
	 tip = dmy \ DMcommand;
	 DMcommand = DMcommand - dmy * tip; 
	 % cmd = randn(97,1) * 0.08; 
	 subplot(1,3,3); 
    bar(DMcommand); 
    title('DM control signals'); 
	 dmctrl.Data = single(DMcommand); 
    
    pause(0.01); 
end
% now, need to convert this to physical units -- microns or waves, 
% based on the curvature of the point source. 