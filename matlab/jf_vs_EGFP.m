% to measure the bleaching of JF-669 at various wavelengths. 
N = 20; 
dat = zeros(N, 2);
wavelength = input('wavelength'); 
power = input('power');
jfconc = input('JFconc');

udpr = dsp.UDPReceiver('LocalIPPort', 31313,'MessageDataType','double'); 

setup(udpr); 

% drain the buffer, if need be. 
datarx = 1; 
while ~isempty(datarx)
	datarx = udpr(); 
end

disp('grabbing data..'); 
tic

for i = 1:N
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
	sumstd
	dat(i, 1) = sumstd; 
	dat(i, 2) = toc; 
end
plot(dat(:,2), dat(:,1))
save(['jf669_vs_EGFP_' num2str(wavelength) '_' num2str(jfconc*10) '.mat'], 'dat', 'wavelength', 'power')