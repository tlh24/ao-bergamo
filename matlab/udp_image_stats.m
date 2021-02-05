close all

global DMcommand dmx dmy dmctrl mmf sock dmangle; 
global mask save_frames save_dmcommand save_wfs_dx save_wfs_dy save_sumstd
global DMcommandHist DMcommandStd DMcommandK
global temperature temperatures k

N = 15e3; 
save_frames = single(zeros(N, 256, 256)); % just to be double sure. 
save_dmcommand = single(zeros(N, 97));
save_sumstd = single(zeros(N, 10)); 
save_wfs_dx = single(zeros(N, sum(mask)));  
save_wfs_dy = single(zeros(N, sum(mask))); 

mmf = memmapfile('../shared_centroids.dat','Format','single','Offset',0,'Repeat',6000);

dmctrl = memmapfile('../shared_dmctrl.dat','Format','single','Offset',0,'Repeat',97, 'Writable',true);

if 1
	DMcommand = zeros(97, 1); 
else
	load('../Best_DMcommand4.mat'); 
	Best_DMcommand = reshape(Best_DMcommand, 97, 1); % transpose
	DMcommand = Best_DMcommand; 
end
DMcommandHist = repmat(DMcommand, 1, 100); 
DMcommandStd = zeros(1, 100);
DMcommandK = zeros(1, 100); 

[dmx, dmy] = dm_actuator_to_xy();
dmangle = angle([dmx + 1i*dmy]); 
dmctrl.Data = single(DMcommand); 

sock = tcpip('0.0.0.0', 31313, 'NetworkRole', 'server');
sock.InputBufferSize = 512*513; 
disp('ok go.'); 
fopen(sock); % waits for a connection 


temperature = 0.01; 
starttemp = 0.0075; % naive start = 0.008
endtemp = 0.000; 
temperatures = linspace(starttemp, endtemp, N); 
k = 1; 

while k < N
	temperature = temperatures(k); 
	pick = floor(rand(1) * 99) + 1; 
	father = DMcommandHist(:,pick); 
	pick2 = pick; 
	while pick2 == pick
		pick2 = floor(rand(1) * 99) + 1; 
	end
	mother = DMcommandHist(:,pick2); 
	recomb = (rand(1)-0.5) * 2*pi; 
	kid = father .* (dmangle > recomb) + mother .* (dmangle <= recomb); 
	noise = (randn(97,1)*temperature) .* (rand(97, 1) > 0.8); 
	send_cmd_get_data(kid+noise)
	
	% this doesn't work. 
	% there are no clean gradients (apparently) in the search space
% 	if mod(k, 1000) == 0
% 		% do a regression
% 		A = save_dmcommand(k-999:k, :); 
% 		A_mean = mean(A, 1); 
% 		A = A - A_mean; 
% 		B = sum(save_sumstd(k-999:k, 2:end), 2); % FIXME!!
% 		B = zscore(B); 
% 		coef = A\B; 
% 		pred = A*coef; 
% 		cc = corrcoef(B, pred)
% 		if cc(1, 2) > 0.2 % noise level is 0.31 about
% 			disp('doing beam search'); 
% 			scoef = coef * std(std(A, [], 2));
% 			% beam search along this direction.
% 			for j = 1:20
% 				cmd = A_mean' + scoef*(j-1)/4; 
% 				send_cmd_get_data(cmd); 
% 			end
% 		end
% 	end
	k = k + 1; 
end

fclose(sock); 

function send_cmd_get_data(cmd)
	global DMcommand dmx dmy DMcommandP dmctrl mmf sock; 
	global mask save_frames save_dmcommand save_wfs_dx save_wfs_dy save_sumstd
	global DMcommandHist DMcommandStd DMcommandK
	global temperature k

	DMcommand = reshape(cmd, 97, 1); 
	DMcommand = min(DMcommand, 0.15); % clip, per reality.
	DMcommand = max(DMcommand, -0.15); 
	DMcommandP = DMcommand - mean(DMcommand); %piston
	tilt = dmx \ DMcommandP;
	% there still could be tip-tilt in the command, 
	% even if it was removed from the WFS derivatives.
	DMcommandP = DMcommandP - dmx * tilt; 
	tip = dmy \ DMcommandP;
	DMcommandP = DMcommandP - dmy * tip; 
	dmctrl.Data = single(DMcommandP); 
	% C++ control program drive the mirrors. 
	
	framen = 0; 
	sumstd = 0; 
	sumframe = single(zeros(256, 256));
	while framen < 10
		datarx = fread(sock, 3, 'double');
		framemin = datarx(2); 
		framemax = datarx(3);
		frame = fread(sock, 256*256, 'uint8'); 
		frame = single(reshape(frame, 256, 256)); 
		frame = (frame / 255.0) * (framemax-framemin) + framemin; 
		save_sumstd(k, framen+1) = datarx(1); 
		if framen > 0
			sumstd = sumstd + datarx(1); 
			sumframe = sumframe + frame; 
		end
		framen = framen+1; 
	end
% 	imagesc(sumframe); 
% 	drawnow; 
	
	% read in the wavefront sensor
	dat = mmf.Data; 
   datx = dat(1:2:end); 
   daty = dat(2:2:end); 
	dx = datx(mask); 
	dy = daty(mask); 
	
	% save some data... 
	save_frames(k, :, :) = single(sumframe); % just to be double sure. 
	save_dmcommand(k, :) = DMcommandP;
	save_wfs_dx(k, :) = dx; 
	save_wfs_dy(k, :) = dy; 
	
	DMcommandHist = [DMcommandHist DMcommandP]; 
	DMcommandStd = [DMcommandStd sumstd]; 
	DMcommandK = [DMcommandK k]; 
	agedecay = 1-((k - DMcommandK) / 5000); 
	[~, indx] = sort(DMcommandStd .* agedecay, 'descend'); 
	DMcommandStd = DMcommandStd(indx(1:100)); 
	DMcommandHist = DMcommandHist(:, indx(1:100)); 
	DMcommandK = DMcommandK(:, indx(1:100)); 
	disp([DMcommandStd(1) mean(DMcommandStd) temperature*1e7]); % display the best one. 
end
