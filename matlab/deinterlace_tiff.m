function [out_stack] = deinterlace_tiff(fname)
% deinterlace tiff based on between scanline correlation. 
% upsamples the image by 4x in the x direction to perform sub-pixel
% alignment. 
tiff_info = imfinfo(fname);

tiff_stack = imread(fname, 1) ; 
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread(fname, ii);
    tiff_stack = cat(3 , tiff_stack, temp_tiff);
end

res = size(tiff_stack, 1);
out_stack = tiff_stack; 
for slice = 1:size(tiff_stack, 3)
	img = tiff_stack(:, :, slice); 
	out = imresize(img, [res res*4]); 
	odd = out(1:2:end, :); 
	even = out(2:2:end, :); 

	for r = 1:res/2
		xc2 = xcorr(even(r,:), odd(r,:));
		if r == 1
			xc = xc2;
		else
			xc = xc + xc2;
		end
		if r < res/2
			xc3 = xcorr(even(r, :), odd(r+1, :)); 
			xc = xc + xc3; 
		end
	end
	[~, col] = max(xc(:));
	offset = col - res*4; 
	% there shouldn't be a one-off error here... 
	% checked with known shift. 

	if offset > 0
		out(2:2:res, 1:end-offset) = even(:, 1+offset:end); 
	end
	if offset < 0
		out(2:2:res, 1-offset:end) = even(:, 1:end+offset);
	end

	out = imresize(out, [res res]); 
	out_stack(:, :, slice) = out; 
	
	if 0
		figure
		subplot(1,2,1)
		imagesc(sqrt(abs(img))); 
		subplot(1,2,2)
		imagesc(sqrt(abs(out))); 
		colormap gray
	end
	disp(['slice ' num2str(slice) ' done']); 
end