#include <pylon/PylonIncludes.h>
#include <math.h>
#include <malloc.h>
#include <png.h> 
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <mutex>
#include <condition_variable>
#include <random>

#include <gtk/gtk.h>
#include <gdk/gdk.h>

#include "gettime.h"
#include <matio.h>
#include "serialize.h"
#include "globals.h"

// Namespace for using pylon objects.
using namespace Pylon;

// Namespace for using GenApi objects.
using namespace GenApi;

// Namespace for using cout.
using namespace std;


// need to take in locations of the flat-wavefront pixel centers
// and convert that to a list of nearest pixels
// which are then summed to calculate the center of mass (perhaps including a threshold... 

float g_centroids[MAX_LENSLETS][2];
float g_centroidsCalib[MAX_LENSLETS][2];
int g_lensletStarts[MAX_LENSLETS][2];
int g_nCentroids = 0;
int g_lensSpanW = 35; 
int g_lensSpanH = 35;
int g_nFrames = 0;
bool g_calibrated = false; 

float g_pixthresh = 70.0/255.0; 
float g_sumthresh = 3.0; 

int g_exposure = 50; 
int g_actuator = 0; 
bool g_set_exposure = true; 
bool g_reset_data = false;
bool g_record_data = false; 
bool g_write_data = false; 

int g_dm_actuator = 0; 
int g_dm_counter = 0;

void init_centroids(const unsigned char* image, 
						  int sx, int sy){
	//fill out the centroid structure in raster-scanning order starting from sx,sy
	int dx = g_lensSpanW; 
	int dy = g_lensSpanH; 
	int n = 0; 
	int xdir = 1; 
	while(sx < g_w-dx/2 && sy < g_h-dy/2 && n < MAX_LENSLETS){
		float cx, cy, cd; 
		cx = cy = cd = 0.0; 
		for(int y = sy - dy/2; y < sy + dy/2; y++){
			for(int x = sx - dx/2; x < sx + dy/2; x++){
				float d = float(image[y*g_w + x]) / 255.f; 
				if(d > g_pixthresh - 0.22f){ // about 0.15 for 15x sigmoid
					float w = (tanhf(( d - g_pixthresh ) * 15.f) + 1.f)/2.f; 
					d *= w; 
					cx += x * d; 
					cy += y * d; 
					cd += d; 
				}
			}
		}
		if(cd > g_sumthresh){
			cx /= cd; 
			cy /= cd; 
			g_lensletStarts[n][0] = cx - dx/2; 
			g_lensletStarts[n][1] = cy - dy/2; 
			printf("centroid %d found at %f %f val %f\n", n, cx, cy, cd); 
			n++; 
			sx = round(cx); 
			sy = round(cy); 
		} else {
			cx /= cd; 
			cy /= cd; 
			printf("rejected centroid at %f %f val %f expected %d %d \n", cx, cy, cd, sx, sy); 
		}
		sx += dx * xdir; 
		if(sx >= g_w-dx/2 || sx < dx/2){
			sx -= dx * xdir; 
			xdir *= -1; 
			sy += dy; 
		}
	}
	printf("found %d lenslets.\n", n); 
	g_nCentroids = n; 
	g_lensSpanW = dx; 
	g_lensSpanH = dy; 
}
void init_centroids_default(const unsigned char* image){
	init_centroids(image, 319, 122);
}

struct calc_centroids_data {
	int	startLens; 
	int	endLens; 
	const unsigned char* image; 
	float* out; 
} ; 

void* calc_centroids(void* v){
	calc_centroids_data * p = (calc_centroids_data*) v; 
	// will perhaps eventually do this on the gpu
	// 'out' needs to be pre-allocated.
	for(int i=p->startLens; i<p->endLens && i<g_nCentroids; i++){
		float cx, cy, cd; 
		cx = cy = cd = 0.0;
		for(int row = 0; row < g_lensSpanH; row++){
			int rr = row + g_lensletStarts[i][1]; 
			for(int col = 0; col < g_lensSpanW; col++){
				int cc = col + g_lensletStarts[i][0]; 
				float d = (float)(p->image[rr * g_w + cc]) / 255.f; 
				if(d > g_pixthresh - 0.22f){ // about 0.15 for 15x sigmoid
					//hard thresholding leads to hard transitions
					//add a secondary sigmoidal weighting around the threshold. 
					float w = (tanhf(( d - g_pixthresh ) * 15.f) + 1.f)/2.f; 
					d *= w; 
					// further weight the COM by pixel intensities. 
					cx += (float)cc * d; 
					cy += (float)rr * d; 
					cd += d; 
				}
			}
		}
		if(cd > g_sumthresh){
			cx /= cd; 
			cy /= cd; 
			p->out[i*2 + 0] = cx; 
			p->out[i*2 + 1] = cy; 
		}
	}
	return NULL; 
}


int writeGrayPNG(const char* filename, int width, int height, const unsigned char *buffer)
{
	int code = 0;
	FILE *fp = NULL;
	png_structp png_ptr = NULL;
	png_infop info_ptr = NULL;
	png_bytep row = NULL;
	
	// Open file for writing (binary mode)
	fp = fopen(filename, "wb");
	if (fp == NULL) {
		fprintf(stderr, "Could not open file %s for writing\n", filename);
		code = 1;
		goto finalise;
	}

	// Initialize write structure
	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if (png_ptr == NULL) {
		fprintf(stderr, "Could not allocate write struct\n");
		code = 1;
		goto finalise;
	}

	// Initialize info structure
	info_ptr = png_create_info_struct(png_ptr);
	if (info_ptr == NULL) {
		fprintf(stderr, "Could not allocate info struct\n");
		code = 1;
		goto finalise;
	}

	// Setup Exception handling
	if (setjmp(png_jmpbuf(png_ptr))) {
		fprintf(stderr, "Error during png creation\n");
		code = 1;
		goto finalise;
	}

	png_init_io(png_ptr, fp);

	// Write header (8 bit colour depth)
	png_set_IHDR(png_ptr, info_ptr, width, height,
			8, PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE,
			PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

	png_write_info(png_ptr, info_ptr);

	// Write image data
	int x, y;
	for (y=0 ; y<height ; y++) {
		png_write_row(png_ptr, &(buffer[y*width]));
	}

	// End write
	png_write_end(png_ptr, NULL);

	finalise:
	if (fp != NULL) fclose(fp);
	if (info_ptr != NULL) png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
	if (png_ptr != NULL) png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
	if (row != NULL) free(row);

	return code;
}

unsigned char* read_png_file(const char *filename) {
	png_byte color_type;
	png_byte bit_depth;
	png_bytep *row_pointers = NULL;
	
	FILE *fp = fopen(filename, "rb");

	png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if(!png) abort();

	png_infop info = png_create_info_struct(png);
	if(!info) abort();

	if(setjmp(png_jmpbuf(png))) abort();

	png_init_io(png, fp);

	png_read_info(png, info);

	g_w      = png_get_image_width(png, info);
	g_h     = png_get_image_height(png, info);
	color_type = png_get_color_type(png, info);
	bit_depth  = png_get_bit_depth(png, info);
	
	printf("reading %s size %d %d depth %d\n", filename, g_w, g_h, bit_depth); 

	// Read any color_type into 8bit depth, RGBA format.
	// See http://www.libpng.org/pub/png/libpng-manual.txt

	if(bit_depth == 16)
		png_set_strip_16(png);

	if(color_type == PNG_COLOR_TYPE_PALETTE)
		png_set_palette_to_rgb(png);

	// PNG_COLOR_TYPE_GRAY_ALPHA is always 8 or 16bit depth.
	if(color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8)
		png_set_expand_gray_1_2_4_to_8(png);

	if(png_get_valid(png, info, PNG_INFO_tRNS))
		png_set_tRNS_to_alpha(png);

	// These color_type don't have an alpha channel then fill it with 0xff.
	if(color_type == PNG_COLOR_TYPE_RGB ||
		color_type == PNG_COLOR_TYPE_GRAY ||
		color_type == PNG_COLOR_TYPE_PALETTE)
		png_set_filler(png, 0xFF, PNG_FILLER_AFTER);

	if(color_type == PNG_COLOR_TYPE_GRAY ||
		color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
		png_set_gray_to_rgb(png);

	png_read_update_info(png, info);

	if (row_pointers) abort();

	row_pointers = (png_bytep*)malloc(sizeof(png_bytep) * g_h);
	for(int y = 0; y < g_h; y++) {
		row_pointers[y] = (png_byte*)malloc(png_get_rowbytes(png,info));
	}

	png_read_image(png, row_pointers);
	png_read_end(png, NULL); 
	fclose(fp);
	
	unsigned char* grayb = (unsigned char*)malloc(g_w * g_h); 
	int n = 0; 
	for(int y=0; y<g_h; y++){
		for(int x=0; x<g_w; x++){
			grayb[n] = ((png_bytep)(row_pointers[y]))[x*4 + 1]; 
			n++; 
		}
	}
	png_destroy_read_struct(&png, &info, NULL);
	
	return grayb; 
}

bool read_calibration_flat(){
	// read in the calibration matlab file. 
	mat_t    *matfp; 
	matvar_t *matvar; 
	bool calibrated = false;

	matfp = Mat_Open("calibration_flat.mat",MAT_ACC_RDONLY); 
	if ( NULL == matfp ) { 
		fprintf(stderr,"Error opening calibration_flat.mat"); 
		return false; 
	} 

	matvar = Mat_VarRead(matfp,"calib"); 
	if ( NULL != matvar ) { 
		if(matvar->rank == 2 && matvar->dims[1] == 2){
			int n = matvar->dims[0]; 
			g_nCentroids = n;
			if(matvar->data_type == MAT_T_DOUBLE){
				for(int i=0; i<MAX_LENSLETS*2; i++)
					g_centroidsCalib[0][i] = -1.0; 
				double* dat = (double*)(matvar->data);
				for(int j=0; j<2; j++){
					for(int i=0; i<n; i++){
						//transpose matlab order to C order
						g_centroidsCalib[i][j] = dat[i + n*j]; 
					}
				}
				printf("read in %d calibrated centroids!\n", n); 
				calibrated = true; 
			} else {
				printf("calibration_flat.m variable calib is not type double\n"); 
			}
		} else {
			printf("calibration_flat.m variable calib needs to be 3000x2\n"); 
		}
		Mat_VarFree(matvar); 
	} 
	Mat_Close(matfp); 
	//now, need to re-initialize the lensletStarts structure
	if(calibrated){
		for(int i=0; i<g_nCentroids; i++){
			g_lensletStarts[i][0] = round(g_centroidsCalib[i][0]) - g_lensSpanW/2; 
			g_lensletStarts[i][1] = round(g_centroidsCalib[i][1]) - g_lensSpanH/2; 
			if(g_lensletStarts[i][0] < 0)
				g_lensletStarts[i][0] = 0; 
			if(g_lensletStarts[i][1] < 0)
				g_lensletStarts[i][1] = 0; 
			if(g_lensletStarts[i][0] > 2048 - g_lensSpanW -1)
				g_lensletStarts[i][0] = 2048 - g_lensSpanW -1; 
			if(g_lensletStarts[i][1] > 2048 - g_lensSpanH -1)
				g_lensletStarts[i][1] = 2048 - g_lensSpanH -1; 
		}
	}
	return calibrated;
} 	

void* video_thread(void*){
	int exitCode = 0;
	g_startTime = gettime(); 
	double lastFrameTime = g_startTime; 
	
	g_calibrated = read_calibration_flat(); 

	PylonInitialize();
	
	const uint8_t *pImageBuffer;
	
	VectorSerialize2<float>* centroidVec; 
	centroidVec = new VectorSerialize2<float>(3000, MAT_C_SINGLE); 
	VectorSerialize<float>* dmVec; 
	dmVec = new VectorSerialize<float>(97, MAT_C_SINGLE); 
	
	//share data with matlab through a mem-mapped file.

	int mmfd = open("shared.dat", O_RDWR, S_IREAD | S_IWRITE);
	if (mmfd < 0) {
		perror("Could not open file for memory mapping");
		exit(1);
	}
	float* mmap_ptr = (float*)mmap(NULL, MAX_LENSLETS*2*4, PROT_READ | PROT_WRITE, MAP_SHARED, mmfd, 0);
	
	if (mmap_ptr == MAP_FAILED) {
		perror("Could not memory map file");
		exit(1);
	}
	close(mmfd); 
	
	//distribute centroid calculation over multiple threads
	pthread_t thread[NTHREADS];
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	calc_centroids_data centroids_thread_data[NTHREADS]; 
	int sl = 0; 
	for(int i=0; i<NTHREADS; i++){
		centroids_thread_data[i].startLens = sl; 
		sl += g_nCentroids / NTHREADS; 
		if(i==NTHREADS-1) sl = g_nCentroids; 
		centroids_thread_data[i].endLens = sl; 
		centroids_thread_data[i].out = &(g_centroids[0][0]); 
	}
	
	//start up the deformable mirror thread. 
	float dm_data[97]; 
	for(int k=0; k<97; k++){
		dm_data[k] = 0.f; 
	}
	pthread_t dm_thread;
	pthread_create( &dm_thread, &attr, dmControl_thread, (void*)dm_data); 
	dm_semaphore = new Semaphore(); 
	
	//init gaussian random number generator, too
	std::default_random_engine generator;
	std::normal_distribution<float> distribution(0.0,0.018);
	
	try{
		CInstantCamera camera( CTlFactory::GetInstance().CreateFirstDevice());

		// Print the model name of the camera.
		cout << "Using device " << camera.GetDeviceInfo().GetModelName() << endl;
		
		// Open the camera.
		camera.Open();

		// Get the required camera settings.
		CIntegerParameter width( camera.GetNodeMap(), "Width");
		CIntegerParameter height( camera.GetNodeMap(), "Height");
		CEnumParameter pixelFormat( camera.GetNodeMap(), "PixelFormat");

		// Map the pixelType
		CPixelTypeMapper pixelTypeMapper(&pixelFormat);
		EPixelType pixelType = pixelTypeMapper.GetPylonPixelTypeFromNodeValue(pixelFormat.GetIntValue());

		// Start the grabbing of c_countOfImagesToGrab images.
		// The camera device is parameterized with a default configuration which
		// sets up free running continuous acquisition.
		camera.StartGrabbing(GrabStrategy_LatestImages, GrabLoop_ProvidedByUser);

		cout << "Please wait. Images are being grabbed." << endl;

		// This smart pointer will receive the grab result data.
		CGrabResultPtr ptrGrabResult;
		// I assume this does not need to be deallocated? 

		while ( camera.IsGrabbing() && !g_die)
		{
			// Wait for an image and then retrieve it. A timeout of 5000 ms is used.
			camera.RetrieveResult( 5000, ptrGrabResult, TimeoutHandling_ThrowException);

			uint8_t *pImageBuffer = (uint8_t *) ptrGrabResult->GetBuffer();

			if(pImageBuffer){
				memcpy(g_data[0], pImageBuffer, g_w*g_h); 
				g_copy[0] = 1; 
			}
			if(!g_calibrated){
				//write out a png of the image!
				writeGrayPNG("framecap.png", 
									ptrGrabResult->GetWidth(), 
									ptrGrabResult->GetHeight(), 
									pImageBuffer); 
				init_centroids_default(pImageBuffer); 
				centroidVec->clear(); 
				g_calibrated = true; 
			} else {
				double start = gettime(); 
				for(int i=0; i<NTHREADS; i++){
					centroids_thread_data[i].image = pImageBuffer; 
					pthread_create( &thread[i], &attr, calc_centroids, (void*)(&(centroids_thread_data[i]))); 
				}
				void* ptr; 
				for(int i=0; i<NTHREADS; i++){
					pthread_join(thread[i], &ptr); 
				}
				g_centroidCalc_label.set(gettime() - start); 
				g_framerate_label.set(1.0 / (start - lastFrameTime));
				g_dataSize_label.set((float)centroidVec->nstored()); 
				lastFrameTime = start; 
				if(g_nFrames%5 == 4){
					if(g_record_data && centroidVec->nstored() < 180000){ 
						for(int i=0; i<g_nCentroids && i<3000; i++){
							centroidVec->m_stor[i] = g_centroids[i][0]; 
							centroidVec->m_stor2[i] = g_centroids[i][1]; 
						}
						for(int i=0; i<97; i++){
							dmVec->m_stor[i] = dm_data[i]; 
						}
						centroidVec->store(); 
						dmVec->store(); 
					}
					for(int i=0; i<97; i++){
						dm_data[i] = 0.f;
					}
					if(0){
						//generate a new dm command signal. 
						for(int i=0; i<97; i++){
							dm_data[i] = distribution(generator);
							if(dm_data[i] > 0.10) dm_data[i] = 0.10; 
							if(dm_data[i] <-0.10) dm_data[i] =-0.10; 
						}
					}
					if(0){
						//need to encode actuator positions into RWC
						float scl = 1.f; 
						if(g_dm_counter & 0x1){
							scl *= -1.f; 
						}
						g_dm_counter++; 
						dm_data[60] = dm_data[49] = dm_data[38] = scl * 0.08; 
						dm_data[59] = dm_data[48] = dm_data[37] = scl * 0.08; 
						dm_data[58] = dm_data[47] = dm_data[36] = scl * 0.08; 
						dm_data[72] = dm_data[61] = dm_data[50] = dm_data[39] = dm_data[28] = scl * 0.05; 
						dm_data[68] = dm_data[57] = dm_data[46] = dm_data[35] = dm_data[24] = scl * 0.05; 
						dm_data[71] = dm_data[70] = dm_data[69] = scl * 0.05; 
						dm_data[27] = dm_data[26] = dm_data[25] = scl * 0.05; 
					}
					if(0){
						if(g_actuator > 96) g_actuator = 96; 
						if(g_actuator < 0) g_actuator = 0; 
						if(g_dm_counter & 0x1){
							dm_data[g_actuator] = -0.10; 
						}else{
							dm_data[g_actuator] = 0.10; 
						}
						g_dm_counter++; 
					}
					dm_semaphore->notify(); 
				}
				memcpy(mmap_ptr, g_centroids, MAX_LENSLETS*2*4); 
			}
			
			if(g_set_exposure){
				INodeMap& nodemap = camera.GetNodeMap();
				// Determine the current exposure time
				double d = CFloatPtr(nodemap.GetNode("ExposureTime"))->GetValue();
				printf("previous exposure time, us: %f\n", d); 
				// Set the exposure time mode to Standard
				// Note: Available on selected camera models only
				CFloatPtr(nodemap.GetNode("ExposureTime"))->SetValue((float)g_exposure);
				g_set_exposure = false; 
			}
			
			if(g_reset_data){
				centroidVec->clear(); 
				g_reset_data = false; 
			}
			
			if(g_write_data){
				vector<Serialize*> g_objs; 
				g_objs.push_back(centroidVec); 
				g_objs.push_back(dmVec); 
				string fname = "centroids.mat";
				writeMatlab(g_objs, (char *)fname.c_str(), false);
				g_write_data = false;
			}
			g_nFrames++; 
		}
	}
	catch (const GenericException &e){
		// Error handling.
		cerr << "An exception occurred." << endl
		<< e.GetDescription() << endl;
		exitCode = 1;
		
		//read in a png and process it. 
		unsigned char* image = read_png_file("framecap.png"); 
		writeGrayPNG("framecap_test.png", g_w, g_h, image); 
		if(image){
			memcpy(g_data[0], image, g_w*g_h); 
			g_copy[0] = 1; 
		}
		long double start = gettime(); 
		init_centroids_default(image); 
		std::cout << "centroid init time: " << gettime() - start << endl; 
		
		start = gettime(); 
		
		for(int i=0; i<NTHREADS; i++){
			centroids_thread_data[i].image = image; 
			pthread_create( &thread[i], &attr, calc_centroids, (void*)(&(centroids_thread_data[i]))); 
		}
		void* ptr; 
		for(int i=0; i<NTHREADS; i++){
			pthread_join(thread[i], &ptr); 
		}
		std::cout << "centroid calc time: " << gettime() - start << endl; 
		//3.5ms with -O3 -- we don't need to move this to the GPU.
		for(int i=0; i<g_nCentroids && i<3000; i++){
			centroidVec->m_stor[i] = g_centroids[i][0]; 
			centroidVec->m_stor2[i] = g_centroids[i][1]; 
		}
		centroidVec->store(); 
		free(image); 
	}
	printf("done.\n"); 
	// Releases all pylon resources. 
	PylonTerminate(); 
	
	free(g_data[0]); 
	g_data[0] = NULL; 
	
	dm_semaphore->notify(); 
	void* ptr; 
	pthread_join(dm_thread, &ptr); 
	
	if (munmap(mmap_ptr, MAX_LENSLETS*8) < 0) {
		perror("Could not unmap file");
		exit(1);
	}
	long exC = (long)exitCode; 
	return (void*)exC;
}
