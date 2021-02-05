%%--- 
load('../calibration_forward.mat'); % valid lenslets in the actual microscope
% load('../calibration_flat.mat');
calib = [mx my]; % this is not actually flat, but just run with it.. 
nlenslets = size(calib, 1); 

mmf = memmapfile('../shared_centroids.dat','Format','single','Offset',0,'Repeat',6000);

dmctrl = memmapfile('../shared_dmctrl.dat','Format','single','Offset',0,'Repeat',97, 'Writable',true);

DMcommand = zeros(97, 1); 
DMcommandSave = DMcommand; 
[dmx, dmy] = dm_actuator_to_xy();
dmctrl.Data = single(DMcommand); 

if(0)
	DMcommandHist = zeros(97, 100); 
	DMcommandErr = ones(1, 100) * 1e10; 
end

fig1 = figure;

ButtonHandle = uicontrol('Style', 'PushButton', ...
                         'String', 'Stop loop', ...
                         'Callback', 'n=10e6');
n=0; 
actuator = 1; 
beamN = 1; 
beamErr = zeros(5, 1); 
beamOffset = [-0.04 0.02 0.0 0.02 0.04]; 
while n < 1e6
    dat = mmf.Data; 
    datx = dat(1:2:end); 
    daty = dat(2:2:end); 
    
    dx = datx(1:nlenslets) - calib(:,1); 
    dy = daty(1:nlenslets) - calib(:,2); 
    
	 if(0)
		% need to recalulate the mask each time, because we may be adjusting
		% system alignment. 
		mask = (abs(dx) < 10) .* (abs(dy) < 10); 
		mask = mask > 0; 
		dx = dx(mask); 
		dy = dy(mask); 

		cx = mean(calib(mask, 1)); 
		cy = mean(calib(mask, 2)); 
		calibx = calib(mask, 1); 
		caliby = calib(mask, 2); 

		calibx_n = calibx - cx; 
		caliby_n = caliby - cy; 
		calibx_n = (calibx_n / 925.0); 
		caliby_n = (caliby_n / 925.0); 
		nlenslets_n = sum(mask); 
		[Z,dZx,dZy] = zernikes_and_derivatives_cartesian_OSA (calibx_n, caliby_n, 7, "NaN");
	 else
		 % fixed alignment, use the calibration file 
		 % variables 'C' and 'mask'. 
		% dx = dx(mask); 
		% dy = dy(mask); 
		if(n < 5)
			cx = mean(calib(:, 1)); 
			cy = mean(calib(:, 2)); 
			calibx = calib(:, 1); 
			caliby = calib(:, 2); 

			calibx_n = calibx - cx; 
			caliby_n = caliby - cy; 
			calibx_n = (calibx_n / 925.0); 
			caliby_n = (caliby_n / 925.0); 
			nlenslets_n = sum(mask); 
			[Z,dZx,dZy] = zernikes_and_derivatives_cartesian_OSA (calibx_n, caliby_n, 7, "NaN");
		end
		dxbad = abs(dx) > 10; 
		dybad = abs(dy) > 10; 
		dx(dxbad) = 0; % this only works so long as we're trying to flatten.
		dy(dybad) = 0; 
	 end
    % use these derivatives to calculate Zernike coefficients. 

    % serially regress dx and dy to the cartesian derivatives. 
    coef = zeros(36, 1); 
    B = [dx ; dy]; 
    indx = 2; 
    len = 1; 
    while(indx < 36)
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
	 coef2(2:3) = 0; %remove tip/tilt, they are a function of alignment.
    D = Z * coef2; 
	 
	 % now, based on dx and dy, calculate updated control signals .. 
    A = [dx - mean(dx); dy - mean(dy); 1]; %remove tip/tilt
    cmd = A' * Cforward ; 
	 cmd = cmd' * -0.08; 
	 if(1)
		DMcommand = 0.9999*DMcommand + cmd; 
		DMcommand = min(DMcommand, 0.15); % clip, per reality.
		DMcommand = max(DMcommand, -0.15); 
	 end
	 if(0)
		 err = sum(abs(D)); 
		 beamErr(beamN) = err; 
		 if(beamN == 5)
			 [~, best] = min(beamErr); 
			 DMcommand(actuator) = DMcommandSave(actuator) + beamOffset(best); 
			 DMcommand = DMcommand - mean(DMcommand); %piston
			 tilt = dmx \ DMcommand;
			 DMcommand = DMcommand - dmx * tilt; 
			 tip = dmy \ DMcommand;
			 DMcommand = DMcommand - dmy * tip; 
			 DMcommandSave = DMcommand; 
			 beamN = 1; 
			 actuator = floor(rand(1) * 96) + 1; 
		 else
			 beamN = beamN + 1;
		 end
		 DMcommand(actuator) = DMcommandSave(actuator) + beamOffset(beamN); 
	 end
	 if(0)
		err = sum(D.^2);
		DMcommandHist = [DMcommandHist DMcommandP]; 
		DMcommandErr = [DMcommandErr err]; 
		[~, indx] = sort(DMcommandErr); 
		DMcommandErr = DMcommandErr(indx(1:100)); 
		DMcommandHist = DMcommandHist(:, indx(1:100)); 
		pick = floor(rand(1) * 99) + 1; 
		DMcommand = DMcommandHist(:,pick) + randn(97, 1) * 0.0015; 
		DMcommand = min(DMcommand, 0.15); % clip, per reality.
		DMcommand = max(DMcommand, -0.15); 
	 end
	 DMcommandP = DMcommand - mean(DMcommand); %piston
	 tilt = dmx \ DMcommandP;
	 % there still could be tip-tilt in the command, 
	 % even if it was removed from the WFS derivatives.
	 DMcommandP = DMcommandP - dmx * tilt; 
	 tip = dmy \ DMcommandP;
	 DMcommandP = DMcommandP - dmy * tip; 
	 dmctrl.Data = single(DMcommandP); 
	 % send the command out before drawing so the system has a chance to
	 % respond.
	 
	 
    figure(fig1); 
    subplot(2,2,1); 
	 % need another way of displaying the wavefront.
    scatter3(calibx_n, -caliby_n, D); 
	 subplot(2,2,2); 
	 colors = jet(200); 
	 colorindx = round(((D - min(D)) / (max(D) - min(D))) * 200); 
	 colorindx = min(colorindx, 200); 
	 colorindx = max(colorindx, 1); 
	 scatter(calibx_n, -caliby_n, 50, colors(colorindx, :), 'filled'); 
    title('Reconstruction')
    subplot(2,2,3); 
    bar(coef); 
    title('Zernike coefficients'); 
	
	 % cmd = randn(97,1) * 0.08; 
	 subplot(2,2,4); 
	 colorindx = round((DMcommandP + 0.15) / 0.3 * 200); 
	 colorindx = min(colorindx, 200); 
	 colorindx = max(colorindx, 1); 
	 scatter(dmx, dmy, 250, colors(colorindx, :), 'filled')
    title('DM control signals'); 
	 drawnow; 
	 n = n+1; 
end