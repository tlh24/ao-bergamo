% udpr = dsp.UDPReceiver('LocalIPPort', 31313,...
% 	'MessageDataType','double'); 
tic; 
N = 18e3; 
NN = N*1;
bleach_correct = 5000; 
nsamp = 40; 

mmf = memmapfile('../shared_centroids.dat','Format','single','Offset',0,'Repeat',6000);
dmctrl = memmapfile('../shared_dmctrl.dat','Format','single','Offset',0,'Repeat',97, 'Writable',true);

load('../data/calibration_forward.mat','cmask','Cforward','V'); 
% save_frames = single(zeros(N, 256, 256)); % just to be double sure. 
save_dmcommand = single(zeros(NN, 97));
save_kid = single(zeros(NN, 40)); 
save_sumstd = single(zeros(NN, 40)); 
save_wfs_dx = single(zeros(NN, sum(cmask)));  
save_wfs_dy = single(zeros(NN, sum(cmask))); 
save_time = single(zeros(NN, nsamp)); 

[dmx, dmy] = dm_actuator_to_xy();
dmangle = angle([dmx + 1i*dmy]); 
DMcommand = zeros(97, 1);
dmctrl.Data = single(DMcommand); 

Vdm = V' * Cforward; 
% vectors are in Vdm rows. 
% only optimize the first 40 components. minus defocus.
Vdm = Vdm(2:41, :); 

temperature = 1; 
starttemp = 4; % start from zero DM command: 0.005
endtemp = 1.0; 
temperatures = linspace(starttemp, endtemp, N); 
k = 1; 

load('../data/calibration_960nm_PSbeads_long4_geneopt.mat', ...
	'Best_DMcommand'); 

sock = tcpip('0.0.0.0', 18080, 'NetworkRole', 'server');
disp('ok go.'); 
fopen(sock); % waits for a connection 
datarx = fread(sock, 3, 'double');

while k < NN
	% periodically reset the algorithm, 
	% to get new draws from the optimization distro. 
	if mod(k, N) == 1
		kid = zeros(40, 1); 
		kidHist = repmat(kid, 1, 100); 
		DMcommandStd = zeros(1, 100);
		DMcommandK = zeros(1, 100); 
		datarx = fread(sock, 3, 'double');
	end
	
	temperature = temperatures(mod(k-1,N)+1); 
	if k < 100
		noiseax = 1; 
	else
		noiseax = mod(floor((k-100)/100), 40) + 1; 
	end
	pick = floor(rand(1) * 99) + 1; 
	father = kidHist(:,pick); 
	pick2 = pick; 
	while pick2 == pick
		pick2 = floor(rand(1) * 99) + 1; 
	end
	mother = kidHist(:,pick2); 
	recomb = round(rand(40,1)); 
	kid = father .* recomb + mother .* (1-recomb); 
	kid(noiseax) = kid(noiseax) + randn(1) * temperature; 
	
	DMcommand = Vdm' * kid + Best_DMcommand'; 

	DMcommand = min(DMcommand, 0.17); % clip, per reality.
	DMcommand = max(DMcommand, -0.17); 
	DMcommandP = DMcommand - mean(DMcommand); %piston
	tilt = dmx \ DMcommandP;
	% there still could be tip-tilt in the command, 
	% even if it was removed from the WFS derivatives.
	DMcommandP = DMcommandP - dmx * tilt; 
	tip = dmy \ DMcommandP;
	DMcommandP = DMcommandP - dmy * tip; 

	dmctrl.Data = single(DMcommandP); 
	% C++ control program (shwfs) drives the mirror. 

	framen = 0; 
	sumstd = 0; 
% 		sumframe = single(zeros(256, 256));
	while framen < nsamp
		datarx = fread(sock, 3, 'double');
		framemin = datarx(2); 
		framemax = datarx(3); 
		save_sumstd(k, framen+1) = datarx(1); 
		save_time(k, framen+1) = toc; 
		if framen > nsamp/2
			sumstd = sumstd + datarx(1); 
% 				sumframe = sumframe + frame; 
		end
		framen = framen+1; 
	end

	% read in the wavefront sensor
	dat = mmf.Data; 
	datx = dat(1:2:end); 
	daty = dat(2:2:end); 
	dx = datx(cmask); 
	dy = daty(cmask); 

	% save some data... 
	% save_frames(k, :, :) = single(sumframe); % just to be double sure. 
	save_dmcommand(k, :) = DMcommandP;
	save_kid(k, :) = kid; 
	save_wfs_dx(k, :) = dx; 
	save_wfs_dy(k, :) = dy; 
	if k < 15
		% ignore the startup transient
		sumstd = 0; 
	end

	kidHist = [kidHist kid]; 
	DMcommandStd = [DMcommandStd sumstd]; 
	DMcommandK = [DMcommandK k]; 
	agedecay = 1-((k - DMcommandK) / bleach_correct); 
	[~, indx] = sort(DMcommandStd .* agedecay, 'descend'); 
	DMcommandStd = DMcommandStd(indx(1:100)); 
	kidHist = kidHist(:, indx(1:100)); 
	DMcommandK = DMcommandK(:, indx(1:100)); 
	disp([DMcommandStd(1) mean(DMcommandStd) temperature*1e5 sumstd noiseax*1e6]); % display the best one. 

	k = k + 1; 
end
	
fname = ['../rundata/DMoptimization_960nm_mouse5.mat']
save(fname, '-v7.3', 'save_*');

fclose(sock); 
