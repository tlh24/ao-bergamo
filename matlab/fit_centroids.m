%%--- 
clear all
close all
load('calibration_flat.mat'); 
load('centroids_50hz_tip-tilt.mat'); 
 
nlenslets = size(calib, 1); 
x = double(x(1:nlenslets, 1:5)); 
y = double(y(1:nlenslets, 1:5)); 
% plot the variance vs. centroid number. 
vr = std(x, 0, 2) + std(y, 0, 2); 
figure; 
plot(vr); 
title(["std of centroid loc, n = " num2str(size(x,2))]); 
xlabel('centroid number'); 
ylabel('std, pixels'); 
% yeah, noisy, esp around the periphery. 

% select 'valid' centroids based on calibration data. 
cx = mean(calib(:, 1)); 
cy = mean(calib(:, 2)); 
pupil = sqrt((calib(:,1) - cx).^2 + (calib(:,2) - cy).^2) < 950; 
calibx = calib(pupil, 1); 
caliby = calib(pupil, 2); 

hold on
vr2 = vr .* pupil; % mask off non-pupil elements. 
plot(vr2, 'r')
all_avg_std = mean(vr); 
pupil_avg_std = sum(vr2) / sum(pupil); 
legend(['all centroids, \sigma=' num2str(all_avg_std)],['pupil centroids, \sigma=' num2str(pupil_avg_std)])

% looks like 0.045 pixels standard deviation.  
% if 10 pixels gives you 95um sag, 
% and we have ~ 0.05 pixels noise, 
% then equivalent sag noise should be 0.47um
% since we're getting ~ 0.4um, extra averaging from 
% Zernike fit must be improving resolution somewhat.  
% To improve resolution further to 0.1um, need 4x improvement, which means
% we'll need to sample 16 uncorrelated measurements.  
% lets check the average autocorrelation function of the sensor. 

acor = [];
for i = 1:nlenslets
    acr = xcorr([zscore(x(i,:))' zscore(y(i,:))']); 
    acor(:,i) = mean(acr, 2); 
end
figure; 
plot(([1:size(acr,1)] - size(acr,1)/2)*0.02, acor(:, pupil)); 
hold on
plot(([1:size(acr,1)] - size(acr,1)/2)*0.02, mean(acor(:, pupil), 2), 'k', 'Linewidth', 5); 
xlabel('time, seconds')
ylabel('autocorrelation, arb. units')
title('autocorrelation of pupil centroid locations'); 
legend('black = mean autocorr'); 

mx = mean(x, 2); 
my = mean(y, 2); 
dx = mx(pupil) - calibx; 
dy = my(pupil) - caliby; 

figure; 
scatter(dx, dy); 
title('scatter of centroid dx dy'); 
% still pretty noisy. 

% use these derivatives to calculate Zernike coefficients. 
calibx_n = calibx - cx; 
caliby_n = caliby - cy; 
calibx_n = (calibx_n / 950.0); 
caliby_n = (caliby_n / 950.0); 
nlenslets_n = sum(pupil); 
[Z,dZx,dZy] = zernikes_and_derivatives_cartesian_OSA (calibx_n, caliby_n, 7, "NaN");
% serially regress dx and dy to the cartesian derivatives. 
coef = zeros(36, 1); 
B = [dx ; dy]; 
disp(['mean dx ' num2str(mean(B(1:nlenslets_n)))]); 
disp(['mean dy ' num2str(mean(B(nlenslets_n+1:end)))]); 
indx = 2; 
len = 1; 
while(indx < 36)
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
    disp(['mean dx ' num2str(mean(B(1:nlenslets_n)))]); 
    disp(['mean dy ' num2str(mean(B(nlenslets_n+1:end)))]); 
end

% looks alright, I guess. 
% try reconstructing based on the Zernike polynomials. 
Z = squeeze(Z); 
D = Z * coef; 
calibration = 94.7e-6 / 2.3; 
% calibrated sag (defocus) of 94.7um corresponding to 2.3 pixels *integrated*
figure; 
subplot(1,3,1); 
scatter3(calibx_n, caliby_n, D * calibration); 
title('Reconstruction')
subplot(1,3,2); 
bar(coef * calibration); 
title('Zernike coefficients'); 
ylabel('rough scaling to meters')
xlabel('coeficient no')
% plot this with defocus, tilt and tip zeroed out 
subplot(1,3,3); 
coef2 = coef; 
coef2(2:5) = 0; 
E = Z * coef2;  
scatter3(calibx_n, caliby_n, E * calibration); 
title('Reconstruction error')


% ----
% want to run this analysis over the whole 
% centroids_50hz_tip-tilt.mat
% so that we can regress dx and dy to tip and tilt, 
% and see if there are systematic errors / nonlinearities. 
clear all; 
load('calibration_flat.mat'); 
load('centroids_50hz_tip-tilt.mat'); 
nlenslets = size(calib, 1); 
x = double(x(1:nlenslets, :)); 
y = double(y(1:nlenslets, :)); 
% select 'valid' centroids based on calibration data. 
cx = mean(calib(:, 1)); 
cy = mean(calib(:, 2)); 
pupil = sqrt((calib(:,1) - cx).^2 + (calib(:,2) - cy).^2) < 950; 
calibx = calib(pupil, 1); 
caliby = calib(pupil, 2); 
% and normalized coordinates ...
calibx_n = calibx - cx; 
caliby_n = caliby - cy; 
calibx_n = (calibx_n / 950.0); 
caliby_n = (caliby_n / 950.0); 

dx = x(pupil,:) - calibx; 
dy = y(pupil,:) - caliby; 

[Z,dZx,dZy] = zernikes_and_derivatives_cartesian_OSA (calibx_n, caliby_n, 7, "NaN");
Z = squeeze(Z); 
dZx = squeeze(dZx); 
dZy = squeeze(dZy); 
nt = size(x, 2); 
tiptilt = zeros(2, floor(nt/5) + 1); 
for i = [1:5:nt]
    % serially regress dx and dy to the cartesian derivatives. 
    B = [mean(dx(:,i:i+4),2) ; mean(dy(:,i:i+4), 2)]; 
    indx = 2; 
    len = 1; 
    % only regress zernike modes 2 and 3 
    while(indx < 4)
        fin = indx + len; 
        A = []; 
        for i = indx : fin
            A = [A [dZx(:,i); dZy(:,i)]];
        end
        C = A\B; 
        pred = A*C; 
        B = B - pred; 
        if indx == 2
            tiptilt(:, floor(i/5)) = C; 
        end
        indx = fin+1;
        len = len+1;
    end
end

