% need to get the photobleaching rate of the 
% quantum dots for a given laser power and zoom. 

close all
fclose(sock); 
clear sock
sock = tcpip('0.0.0.0', 31313, 'NetworkRole', 'server');
sock.InputBufferSize = 512*513; 
fopen(sock); % waits for a connection 

f1 = figure; 

time = zeros(20000, 1); 
intens = zeros(20000, 1); 
k = 1; 
while k < 20000 
	datarx = fread(sock, 3, 'double');
	frame = fread(sock, 512*512, 'uint8'); 
	frame = reshape(frame, 512, 512); 
	time(k) = toc; 
	intens(k) = datarx(1); 
	disp(['dataReceived: ' num2str(datarx(1)) ' ' num2str(datarx(2))]); 
% 	imagesc(frame); 
% 	drawnow; 
	k = k + 1; 
end

plot(time, intens); 