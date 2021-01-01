%%--- 

load('calibration_flat.mat'); 

nlenslets = size(calib, 1); 
% select 'valid' centroids based on calibration data. 
cx = mean(calib(:, 1)); 
cy = mean(calib(:, 2)); 
pupil = sqrt((calib(:,1) - cx).^2 + (calib(:,2) - cy).^2) < 950; 
calibx = calib(pupil, 1); 
caliby = calib(pupil, 2); 

calibx_n = calibx - cx; 
caliby_n = caliby - cy; 
calibx_n = (calibx_n / 950.0); 
caliby_n = (caliby_n / 950.0); 
nlenslets_n = sum(pupil); 
[Z,dZx,dZy] = zernikes_and_derivatives_cartesian_OSA (calibx_n, caliby_n, 6, "NaN");

mmf = memmapfile('shared.dat','Format','single','Offset',0,'Repeat',6000);

fig1 = figure;

for x = 1:5000
    dat = mmf.Data; 
    dx = dat(1:2:end); 
    dy = dat(2:2:end); 
    
    dx = dx(1:nlenslets) - calib(:,1); 
    dy = dy(1:nlenslets) - calib(:,2); 
    
    dx = dx(pupil); 
    dy = dy(pupil); 

    % use these derivatives to calculate Zernike coefficients. 

    % serially regress dx and dy to the cartesian derivatives. 
    coef = zeros(28, 1); 
    B = [dx ; dy]; 
    indx = 2; 
    len = 1; 
    while(indx < 28)
        fin = indx + len; 
        disp(['working on ' num2str(indx) ' to ' num2str(fin)]); 
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
    D = Z * coef; 
    figure(fig1); 
    subplot(1,2,1); 
    scatter3(calibx_n, caliby_n, D/max(D)); 
    title('Reconstruction')
    subplot(1,2,2); 
    bar(coef); 
    title('Zernike coefficients'); 
    pause(0.01); 
end
% now, need to convert this to physical units -- microns or waves, 
% based on the curvature of the point source. 