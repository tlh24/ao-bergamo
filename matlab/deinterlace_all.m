function [] = deinterlace_all(directory, matlabformat)

olddir = pwd(); 
cd(directory)

a = dir(); 
for j = 1:length(a)
	fname = a(j).name; 
	indx = regexp(fname, 'SUM_'); 
	if numel(indx) > 0
		fname
		out_stack = deinterlace_tiff(fname); 
		if matlabformat
			fname_di = [fname(1:end-4) '_deinterlace.mat']; 
			save(fname_di, 'out_stack')
		else
			fname_di = [fname(1:end-4) '_deinterlace.tif']; 
			t = Tiff(fname_di,'w');
			for i = 1:size(out_stack, 3)
				tagstruct.ImageLength     = size(out_stack,1);
				tagstruct.ImageWidth      = size(out_stack,2);
				tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
				tagstruct.BitsPerSample   = 32;
				tagstruct.SamplesPerPixel = 1;
				tagstruct.RowsPerStrip    = 16;
				tagstruct.SampleFormat    = Tiff.SampleFormat.IEEEFP;
				tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
				tagstruct.Software        = 'MATLAB';
				t.setTag(tagstruct);
				t.write(out_stack(:,:,i))
				t.writeDirectory(); 
			end
			t.close()
		end
	end
end

cd(olddir)
