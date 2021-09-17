% udpr = dsp.UDPReceiver('LocalIPPort', 31313,...
% 	'MessageDataType','double'); 
tic; 
N = 10e3; 
NN = N*10;
bleach_correct = 7000; 

mmf = memmapfile('../shared_centroids.dat','Format','single','Offset',0,'Repeat',6000);
dmctrl = memmapfile('../shared_dmctrl.dat','Format','single','Offset',0,'Repeat',97, 'Writable',true);

load('../data/calibration_forward.mat', 'cmask'); 
% save_frames = single(zeros(N, 256, 256)); % just to be double sure. 
save_dmcommand = single(zeros(NN, 97));
save_sumstd = single(zeros(NN, 12)); 
save_wfs_dx = single(zeros(NN, sum(cmask)));  
save_wfs_dy = single(zeros(NN, sum(cmask))); 
save_time = single(zeros(NN, 12)); 

[dmx, dmy] = dm_actuator_to_xy();
dmangle = angle([dmx + 1i*dmy]); 
DMcommand = zeros(97, 1);
dmctrl.Data = single(DMcommand); 

% % drain the udp port. 
% datarx = 1; 
% while ~isempty(datarx)
% 	datarx = udpr(); 
% 	disp(['drain buffer size ' num2str(numel(datarx))]); 
% end

temperature = 0.01; 
starttemp = 0.007; % start from zero DM command: 0.005
endtemp = 0.0015; 
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
		if 0
			DMcommand = Best_DMcommand'; 
		else
			DMcommand = zeros(97, 1); 
		end
		DMcommandHist = repmat(DMcommand, 1, 100); 
		DMcommandStd = zeros(1, 100);
		DMcommandK = zeros(1, 100); 
		datarx = fread(sock, 3, 'double');
	end
	
	temperature = temperatures(mod(k-1,N)+1); 
	pick = floor(rand(1) * 99) + 1; 
	father = DMcommandHist(:,pick); 
	pick2 = pick; 
	while pick2 == pick
		pick2 = floor(rand(1) * 99) + 1; 
	end
	mother = DMcommandHist(:,pick2); 
	recomb = (rand(1)-0.5) * 2*pi; 
	kid = father .* (dmangle > recomb) + mother .* (dmangle <= recomb); 
	noise = (randn(97,1)*temperature) .* (rand(97, 1) > 0.83); 

	if 0
		% debug timing: why is the optimization no longer working? 
		kid = zeros(97, 1); 
		if mod(k, 2) == 1 
			noise = zeros(97, 1); 
		else
			noise = randn(97, 1)*temperature; 
		end
	end

	DMcommand = reshape(kid+noise, 97, 1); 
	DMcommand = min(DMcommand, 0.15); % clip, per reality.
	DMcommand = max(DMcommand, -0.15); 
	DMcommandP = DMcommand - mean(DMcommand); %piston
	tilt = dmx \ DMcommandP;
	% there still could be tip-tilt in the command, 
	% even if it was removed from the WFS derivatives.
	DMcommandP = DMcommandP - dmx * tilt; 
	tip = dmy \ DMcommandP;
	DMcommandP = DMcommandP - dmy * tip; 

	% drain the udp port --force sync.
% 	datarx = 1; 
% 	while ~isempty(datarx)
% 		datarx = udpr(); 
% 	end
	dmctrl.Data = single(DMcommandP); 
	% C++ control program drives the mirror. 

	framen = 0; 
	sumstd = 0; 
% 		sumframe = single(zeros(256, 256));
	while framen < 12
		datarx = fread(sock, 3, 'double');
% 		datarx = []; 
% 		while numel(datarx) < 1
% 			datarx = udpr(); 
% 		end
		framemin = datarx(2); 
		framemax = datarx(3);
% 			frame = fread(sock, 256*256, 'uint8'); 
% 			frame = single(reshape(frame, 256, 256)); 
% 			frame = (frame / 255.0) * (framemax-framemin) + framemin; 
		save_sumstd(k, framen+1) = datarx(1); 
		save_time(k, framen+1) = toc; 
		if framen > 6
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
	save_wfs_dx(k, :) = dx; 
	save_wfs_dy(k, :) = dy; 
	if k < 15
		% ignore the startup transitent
		sumstd = 0; 
	end

	DMcommandHist = [DMcommandHist DMcommandP]; 
	DMcommandStd = [DMcommandStd sumstd]; 
	DMcommandK = [DMcommandK k]; 
	if 0
		if mean(k-DMcommandK) > 450
			bleach_correct = bleach_correct * 0.99461 
		end
		if mean(k-DMcommandK) < 110 && k > 1000
			bleach_correct = bleach_correct * 1.00237
		end
	end
	agedecay = 1-((k - DMcommandK) / bleach_correct); 
	% this, roughly, should mirror the photobleaching rate
	% for the given excitation wavelength.
	% seems with 7% power at 950nm, it halves over 10k frames. 
	% or goes from ~ 14 to 10 over 5k frames, slope of ~= 
	[~, indx] = sort(DMcommandStd .* agedecay, 'descend'); 
	DMcommandStd = DMcommandStd(indx(1:100)); 
	DMcommandHist = DMcommandHist(:, indx(1:100)); 
	DMcommandK = DMcommandK(:, indx(1:100)); 
	disp([DMcommandStd(1) mean(DMcommandStd) temperature*1e7 sumstd]); % display the best one. 

	k = k + 1; 
end
	
fname = ['../rundata/DMoptimization_1080nm_orangePSbeads_long1.mat']
save(fname, '-v7.3', 'save_*');

fclose(sock); 
