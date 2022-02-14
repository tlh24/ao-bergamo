[dmx, dmy] = dm_actuator_to_xy();

m = max(sqrt(dmx.^2 + dmy.^2)) + 0.7; 

dmx = dmx / m; 
dmy = dmy / m; 

[Z,dZx,dZy] = zernikes_and_derivatives_cartesian_OSA(...
	dmx, dmy, 7, "NaN");

Z = squeeze(Z); 
dZx = squeeze(dZx); 
dZy = squeeze(dZy);

scl = 0.014; 
scaling = zeros(36, 1); 
% 1 2 4 7 11 16 22
k = 4; 
for order = 3:8
	for mode = 0:order-1
		scaling(k+mode) = scl; 
	end
	k = k + order; 
	scl = scl / 2; 
end

scaling = repmat(scaling, 1, 97); 
Z = Z .* scaling'; 

for k = 1:20
	DMcommandP = Z * randn(36, 1) + randn(97, 1) .* 0.025; 

	subplot(1,2,1)
	colors = jet(200); 
	colorindx = round((DMcommandP + 0.15) / 0.3 * 200); 
	colorindx = min(colorindx, 200); 
	colorindx = max(colorindx, 1); 
	scatter(dmx, dmy, 250, colors(colorindx, :), 'filled')
	title('DM control signals'); 

	subplot(1,2,2)
	stem(DMcommandP); 
	drawnow; 
	pause(0.5); 
end

save('../data/dm_zernike_ctrl.mat', 'Z'); 