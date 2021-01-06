% ----
% want to run this analysis over the whole 
% centroids_50hz_tip-tilt.mat
% so that we can regress dx and dy to tip and tilt, 
% and see if there are systematic errors / nonlinearities. 
clear all; 
close all;
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

if 0
    % lowpass filter the dx and dy signals, see if that makes any
    % difference
    [b,a] = butter(3,0.2,'low');     
    for i = 1:size(dx, 1)
        dxf = filtfilt(b,a,dx(i,:)); 
        dyf = filtfilt(b,a,dy(i,:)); 
        dx(i,:) = dxf; 
        dy(i,:) = dyf; 
    end
end

[Z,dZx,dZy] = zernikes_and_derivatives_cartesian_OSA (calibx_n, caliby_n, 7, "NaN");
Z = squeeze(Z); 
dZx = squeeze(dZx); 
dZy = squeeze(dZy); 
nt = size(x, 2); 
coef = zeros(36, nt); 
rmserr = zeros(8, nt); 
residualsx = zeros(sum(pupil), nt); 
residualsy = zeros(sum(pupil), nt); 
for i = [1:nt]
    % serially regress dx and dy to the cartesian derivatives. 
    %first mask off data that didn't change from last frame.
    if i > 1
        mask = (abs(dx(:,i-1) - dx(:,i)) + abs(dy(:,i-1) - dy(:,i))) > 0.0001; 
    else
        mask = ones(size(dx,1), 1); 
    end
    B = [dx(mask,i) ; dy(mask,i)]; 
    indx = 2; 
    len = 1; 
    cnt = 2; 
    ccoef = zeros(36, 1); 
    rrmserr = zeros(8, 1); 
    rrmserr(1) = std(B); 
    % only regress zernike modes 2 and 3 
    while(indx < 36)
        fin = indx + len; 
        % disp(['working on ' num2str(indx) ' to ' num2str(fin)]);
        A = []; 
        for j = indx : fin
            A = [A [dZx(mask,j); dZy(mask,j)]];
        end
        C = A\B; 
        pred = A*C; 
        B = B - pred; 
        ccoef(indx:fin) = C; 
        rrmserr(cnt) = std(B); 
        indx = fin+1;
        len = len+1;
        cnt = cnt+1; 
    end
    coef(:, i) = ccoef; 
    rmserr(:, i) = rrmserr; 
    bx = B(1:sum(mask)); 
    by = B(sum(mask)+1:end); 
    unmask = find(mask > 0); 
    residualsx(unmask, i) = bx; 
    residualsy(unmask, i) = by; 
end

tiptilt = coef(2:3, :); 

colors = ['r', 'g', 'b', 'c', 'm', 'k']; 
for j = [1:6]
    figure; 
    indx = floor(rand(6,1) * sum(pupil))+1; 
    for k = 1:6
        subplot(2,2,1); 
        scatter(tiptilt(2,:), dx(indx(k), :), colors(k)); 
        hold on
        subplot(2,2,2); 
        scatter(tiptilt(1,:), dy(indx(k), :), colors(k)); 
        hold on
        subplot(2,2,3); 
        plot(dx(indx(k), :), colors(k)); 
        hold on
        subplot(2,2,4); 
        plot(dy(indx(k), :), colors(k)); 
        hold on
    end
end

% it seems that some of the data is garbage -- it wasn't updated by the
% centroid algorithm, therefore should be thrown out.  
ddx = diff(dx, 1, 2); 
ddy = diff(dy, 1, 2); 
% imagesc(abs(ddx) < 0.0001); 

figure; 
indx = floor(rand(60,1) * sum(pupil)) + 1; 
plot(dx(indx, :)');  

calibration = 94.7e-6 / 2.3; 
coef = coef * calibration; 

figure; 
indx = 2; 
len = 1; 
cnt = 1; 
while(indx < 36)
    subplot(3,3,cnt);
    fin = indx + len; 
    plot(coef(indx:fin, :)')
    title(['Zernike coef ' num2str(indx) ' to ' num2str(fin) ]); 
    cnt = cnt + 1; 
    indx = fin+1; 
    len = len + 1; 
end

% more exploratory stuff -- try a weiner filter to predict tip-tilt. 
A = [dx ; dy]'; 
B = tiptilt'; 
Ca = A\B; 
tiptiltw = A*Ca; 
plot(tiptilt'); 
hold on
plot(tiptiltw); 
% this works, but need more data!
% also much easier to do a weiner filter than multiple regressions
% could potentially map from raw centroids to corrected centroids 
% using a set of filters ... 
% that said, the right solution is to get a microlens array 
% with a longer focal length! 