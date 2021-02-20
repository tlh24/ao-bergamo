% perform sequential beam-search along the zernike modes 
% to improve image brightness
close all

zernikectrl = memmapfile('../shared_zernike.dat', ...
	'Format','single','Offset',0,'Repeat',36,'Writable',true); 

zcmd = zeros(36, 1); 
zcmdhist = zeros(36, 500);
zcmdv = zeros(1, 500); 
zernikectrl.Data = single(zcmd); 
%piston, tilt and tip are ignored.

N = 200;
save_zcmd = single(zeros(N, 36));
save_sumstd = single(zeros(N, 10)); 

udpr = dsp.UDPReceiver('LocalIPPort', 31313,'MessageDataType','double'); 

setup(udpr); 

% drain the buffer, if need be. 
datarx = 1; 
while ~isempty(datarx)
	datarx = udpr(); 
end

disp('waiting for warm-up frames.'); 
% the first forty frames are warm-up junk too
while framen < 100
	datarx = []; 
	while(length(datarx) < 1)
		datarx = udpr(); 
	end
	framen = framen+1; 
end

disp('grabbing data..'); 

% some baseline activity here ...
for j = 1:5
	framen = 0; 
	sumstd = 0; 
	while framen < 10
		datarx = []; 
		while(length(datarx) < 1)
			datarx = udpr(); 
		end
		save_sumstd(j, framen+1) = datarx(1);
		if framen > 2
			sumstd = sumstd + datarx(1); 
		end
		framen = framen+1; 
	end
	zcmdhist = [zcmdhist zcmd]; 
	zcmdv = [zcmdv sumstd]; 
	[~, indx] = sort(zcmdv, 'descend'); 
	zcmdhist = zcmdhist(:, indx(1:200)); 
	zcmdv = zcmdv(indx(1:200)); 
end

scl = 0.015;
k = 4; 
j = 1; 
for order = 3:8
	for mode = 0:order-1
		for b = -6:6
			zcmd(k+mode) = b * scl; 
			zernikectrl.Data = single(zcmd); 
			if b == -5
				% pause a bit at the beginning of the beams
				% to let the DM settle
				while framen < 5
					datarx = []; 
					while(length(datarx) < 1)
						datarx = udpr(); 
					end
					framen = framen+1; 
				end
			end
			% now need to read in the intensity.. 
			% it may take a few frames for the DM to respond, 
			% plus transmission latencies etc. 
			% scrap the first 2 reeived intensities
			framen = 0; 
			sumstd = 0; 
			while framen < 10
				datarx = []; 
				while(length(datarx) < 1)
					datarx = udpr(); 
				end
				save_sumstd(j, framen+1) = datarx(1);
				if framen > 1
					sumstd = sumstd + datarx(1); 
				end
				framen = framen+1; 
			end
			save_zcmd(j, :) = zcmd;
			
			zcmdhist = [zcmdhist zcmd]; 
			zcmdv = [zcmdv sumstd]; 
			[~, indx] = sort(zcmdv, 'descend'); 
			zcmdhist = zcmdhist(:, indx(1:200)); 
			zcmdv = zcmdv(indx(1:200)); 
			if sum(zcmd == zcmdhist(:,1))==36 && sumstd == zcmdv(1)
				disp([order mode b ...
					zcmdv(1)/1e6, mean(zcmdv)/1e6])
			end
			zcmd = zcmdhist(:,1); % top so far.. 
			zernikectrl.Data = single(zcmd); 
			j = j+1; 
		end
	end
	k = k + order; 
	scl = scl / 1.75; 
end
disp('done.') 

for i = 1:10
	zernikectrl.Data = single(zeros(size(zcmd)));
	pause(1); 
	zernikectrl.Data = single(zcmd); 
	pause(1); 
end

release(udpr); 