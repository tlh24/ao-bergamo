% perform sequential beam-search along the zernike modes 
% to improve image brightness
close all

zernikectrl = memmapfile('../shared_zernike.dat', ...
	'Format','single','Offset',0,'Repeat',36,'Writable',true); 

zcmd = zeros(36, 1); 
zernikectrl.Data = single(zcmd); 
%piston, tilt and tip are ignored.

udpr = dsp.UDPReceiver('LocalIPPort', 31313,'MessageDataType','double'); 

setup(udpr); 

% drain the buffer, if need be. 
datarx = 1; 
while ~isempty(datarx) 
	datarx = udpr(); 
end

disp('grabbing data..'); 

% some baseline measurements to get the noise level
baseline = zeros(20, 1); 
for j = 1:20
	framen = 0; 
	sumstd = 0; 
	while framen < 10
		datarx = []; 
		while(length(datarx) < 1)
			datarx = udpr(); 
		end
		if framen > 1
			sumstd = sumstd + datarx(1); 
		end
		framen = framen+1; 
	end
	baseline(j) = sumstd; 
end

orders = [3 4 5]; % one-indexed. 
scls = [1 (1/1.75) (1/3.0625)] * 0.015; 
indxs = [4 7 11];
for o = 1:length(orders)
	order = orders(o); 
	scl = scls(o); 
	k = indxs(o); 
	for mode = 0:order-1
		if ((order ~= 3) || (mode ~= 1))
			% ignore spherical, it just refocuses. 
			save_intens = zeros(2, 17); 
			save_cmd = zeros(2, 17, 36); 
			for pass = 1:2
				for b = -8:8
					b2 = b * ((pass-1)*2-1); 
					zcmd(k+mode) = b2 * scl; 
					zernikectrl.Data = single(zcmd); 
					if abs(b2) == 6
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
					% scrap the first 2 received intensities
					framen = 0; 
					sumstd = 0; 
					while framen < 10
						datarx = []; 
						while(length(datarx) < 1)
							datarx = udpr(); 
						end
						if framen > 1
							sumstd = sumstd + datarx(1); 
						end
						framen = framen+1; 
					end
					save_intens(pass, b2+9) = sumstd; 
					save_cmd(pass, b2+9, :) = zcmd; 
				end
			end
			intens2 = sum(save_intens, 1); 
			[mx, indx] = sort(intens2, 'descend'); 
			hold off; 
			plot(intens2); 
			hold on; 
			plot(save_intens(1, :)*2, 'r')
			plot(save_intens(2, :)*2, 'r')
			title(['order ' num2str(order) ' mode ' num2str(mode)]); 
			stb2 = 2.4*std(baseline);
			plot(intens2 + stb2, 'g'); 
			plot(intens2 - stb2, 'g'); 
			drawnow;
			if mx(1) - median(intens2) > stb2
				disp([order mode indx(1)-9 mx(1)/1e6]); 
				zcmd = save_cmd(1, indx(1), :); 
			end
		else
			disp('skipping spherical')
		end
	end
end

disp('done evaluating, now blinking') 
for i = 1:10
	zernikectrl.Data = single(zeros(size(zcmd)));
	pause(1); 
	zernikectrl.Data = single(zcmd); 
	pause(1); 
end

release(udpr); 
disp('compeleted'); 