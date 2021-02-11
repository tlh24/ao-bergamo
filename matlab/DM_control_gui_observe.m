%%--- 
load('../data/calibration_flat.mat');
load('../data/calibration_forward.mat'); % valid lenslets in the actual microscope
load('../data/calibration_geneopt.mat');
nlenslets = sum(mask);
% this is confusing, i know. 
% 'mask' from calibration_flat, converts from all possible 3k centroids to
% those physically present on the MLA. 
% mask corresponds to indexes that shwfs will actually output -- 
% shwfs reads in calibration_flat.mat.
% 'cmask' takes those illuminated on the microscope.  
% as of Feb 4, sum(mask) is 1052; 
% sum(cmask) is 539. 

mmf = memmapfile('../shared_centroids.dat',...
	'Format','single','Offset',0,'Repeat',6000);

dmctrl = memmapfile('../shared_dmctrl.dat',...
	'Format','single','Offset',0,'Repeat',97, 'Writable',true);

zernikectrl = memmapfile('../shared_zernike.dat', ...
	'Format','single','Offset',0,'Repeat',36,'Writable',true); 

DMcommand = single(dmctrl.Data); 
[dmx, dmy] = dm_actuator_to_xy();
Zcommand = zeros(36, 1); 
zernikectrl.Data = single(Zcommand); 

fig1 = figure;

ButtonHandle = uicontrol('Style', 'PushButton', ...
                         'String', 'Stop loop', ...
                         'Callback', 'n=10e6');
n=0; 
while n < 1e6
    dat = mmf.Data; 
    datx = dat(1:2:end); 
    daty = dat(2:2:end); 
    
    dx = datx(cmask) - genecalib(:,1); 
    dy = daty(cmask) - genecalib(:,2); 
    
	 % fixed alignment, use the calibration file 
	 % variables 'C' and 'mask'. 
	if(n < 5)
		cx = mean(genecalib(:, 1)); 
		cy = mean(genecalib(:, 2)); 
		calibx = genecalib(:, 1); 
		caliby = genecalib(:, 2); 

		calibx_n = calibx - cx; 
		caliby_n = caliby - cy; 
		calibx_n = (calibx_n / 715.0); 
		caliby_n = (caliby_n / 715.0); 
		nlenslets_n = sum(mask); 
		[Z,dZx,dZy] = zernikes_and_derivatives_cartesian_OSA (calibx_n, caliby_n, 7, "NaN");
		Z = squeeze(Z); 
		dZx = squeeze(dZx); 
		dZy = squeeze(dZy); 
	end
	dxbad = abs(dx) > 30; 
	dybad = abs(dy) > 30; 
	dx(dxbad) = 0; % this only works so long as we're trying to flatten.
	dy(dybad) = 0; 
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
            A = [A [dZx(:,i); dZy(:,i)]];
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
	 coef2 = coef; 
	 coef2(2:3) = 0; %remove tip/tilt, they are a function of alignment.
    D = Z * coef2; 
	 
	 Zcommand(4) = 0.15 * sin(n / 11); 
	 Zcommand(5) = 0.15 * sin(n / 7); 
	 Zcommand(6) = 0.15 * sin(n / 13);
	 
	 Zcommand(7) = 0.05 * sin(n / 17);
	 Zcommand(8) = 0.05 * sin(n / 19);
	 Zcommand(9) = 0.05 * sin(n / 23);
	 Zcommand(10) = 0.05 * sin(n / 29);
	 
	 Zcommand(11) = 0.02 * sin(n / 31);
	 Zcommand(12) = 0.02 * sin(n / 37);
	 Zcommand(13) = 0.02 * sin(n / 41);
	 Zcommand(14) = 0.02 * sin(n / 43);
	 Zcommand(15) = 0.02 * sin(n / 47);
	 zernikectrl.Data = single(Zcommand); 
	 
	 DMcommandP = dmctrl.Data; 
	 
    figure(fig1); 
    subplot(2,2,1); 
	 % need another way of displaying the wavefront.
    scatter3(calibx_n, -caliby_n, D); 
	 axis([-1 1 -1 1 -1 1]); 
	 subplot(2,2,2); 
	 colors = jet(200); 
	 colorindx = round(((D - min(D)) / (max(D) - min(D))) * 200); 
	 colorindx = min(colorindx, 200); 
	 colorindx = max(colorindx, 1); 
	 scatter(calibx_n, -caliby_n, 50, colors(colorindx, :), 'filled'); 
    title('Reconstruction')
    subplot(2,2,3); 
    bar([coef Zcommand]); 
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