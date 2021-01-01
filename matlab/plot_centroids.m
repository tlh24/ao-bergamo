load centroids.mat
x = x(:, 2:end); % first row is garbage for some reason.
y = y(:, 2:end); 
nlenslets = sum(sum(x, 2) > 1); 
x = x(1:nlenslets, :);
y = y(1:nlenslets, :);
ii = floor(rand(1,5) * nlenslets); 
cols = ['r' 'g' 'b' 'c' 'm'];
for i=1:5
    xm = mean(x(ii(i), :)); 
    ym = mean(y(ii(i), :)); 
    plot(x(ii(i), :)-xm, y(ii(i), :)-ym, cols(i), 'LineWidth', 2); 
    hold on
end

% plot the average movement of all the pixels. 
xm = mean(x, 2); 
ym = mean(y, 2); 

xx = x - xm; 
yy = y - ym; 

xxm = mean(xx, 1); 
yym = mean(yy, 1); 

plot(xxm, yym, 'k', 'LineWidth', 3); 
% individual centroids look like crap -- 
% but this seems noise-free and sensitive. 

% alright, need to make a pupil, normalize the coordinates, and fit the
% derivatives. 

