// closed-loop control of the deformable mirror. 
// see DM_control_minimal.m
#include <matio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <stdio.h>

#include <string>
#include <unistd.h>
#include <mutex>
#include <condition_variable>
#include <random>
#include <gtk/gtk.h>
#include <gdk/gdk.h>
#include "globals.h"

#define geneopt_num 6

const char* geneopt_fnames[] = {
	"data/calibration_950nm_1_geneopt.mat",
	"data/calibration_950nm_2_geneopt.mat", 
	"data/calibration_960nm_PSbeads_1_geneopt.mat", 
	"data/calibration_960nm_Retrobeads_1_geneopt.mat", 
	"data/calibration_1020nm_1_geneopt.mat", 
	"data/calibration_1020nm_2_geneopt.mat", 
};
const char* geneopt_names[] = {
	"950nm #1",
	"950nm #2",
	"960nm PSbeads", 
	"960nm Retrobeads", 
	"1020nm #1",
	"1020nm #2",
};

unsigned char g_cmask[3000]; 
int g_activeCentroids = 0; 
int g_nZernike = 0; 
gsl_matrix* g_cforward = NULL; //forward transform from centroids to dm
gsl_matrix* g_genecalib[geneopt_num]; //'system flat' as determined via genetic algos
gsl_matrix* g_dat = NULL; // copy matrix for condensed centroid data
gsl_matrix* g_zcoef = NULL; //commanded zernike coefficients.
gsl_matrix* g_Z = NULL; 
gsl_matrix* g_dZx = NULL; 
gsl_matrix* g_dZy = NULL; 
gsl_matrix* g_dmcommand = NULL; 
gsl_matrix* g_dmZ = NULL; // Zernike matrix for converting gaussian random numbers to a DM command. see zernike_for_actuators.m

bool read_matfile_variable(
	mat_t* matfp, 
	const char* fname, 
	const char* vname, 
	int rows, int cols, 
	gsl_matrix** mat)
{
	matvar_t *matvar = NULL; 
	matvar = Mat_VarRead(matfp,vname); 
	if (matvar != NULL) {
		if(matvar->rank == 2 && matvar->dims[0] == rows && matvar->dims[1] == cols){
			if(matvar->data_type == MAT_T_DOUBLE){
				if(*mat != NULL)
					gsl_matrix_free(*mat); 
				*mat = gsl_matrix_alloc(rows, cols); 
				double* dat = (double*)(matvar->data);
				//transpose to C-order.
				for(int c=0; c<cols; c++){
					for(int r=0; r<rows; r++){
						gsl_matrix_set(*mat, r, c, dat[c*rows + r]); 
					}
				}
			} else {
				printf("variable %s in %s should be type double\n", vname, fname); 
				return false;
			}
		} else {
			printf("variable %s in %s should be %d x %d \n", vname, fname, rows, cols);
			return false;
		}
		Mat_VarFree(matvar); 
	} else {
		printf("could not load variable %s from %s \n", vname, fname);
		return false; 
	}
	return true; 
}

bool dm_control_init()
{
	// initialize these in the most fault-tolerant order.. 
	if(g_dmcommand) gsl_matrix_free(g_dmcommand); 
	g_dmcommand = gsl_matrix_alloc(97, 1); 
	for(int i=0; i<97; i++){
		gsl_matrix_set(g_dmcommand, i, 0, 0.0); 
	}
	g_nZernike = 36; //hardcoded, bad, iknow iknow.
	mat_t *matfp = Mat_Open("data/dm_zernike_ctrl.mat",MAT_ACC_RDONLY); 
	if ( NULL == matfp ) { 
		fprintf(stderr,"Error opening dm_zernike_ctrl.mat"); 
		return false; 
	}
	bool res = read_matfile_variable(matfp, "data/dm_zernike_ctrl.mat", "Z", 97, g_nZernike, &g_dmZ); 
	if(!res){ Mat_Close(matfp); return res; }
	Mat_Close(matfp);
	
	// note: g_nlenslets should be set after reading in calibration_flat.m. 
	matfp = Mat_Open("data/calibration_forward.mat", MAT_ACC_RDONLY); 
	if ( matfp == NULL ) { 
		fprintf(stderr,"Error opening calibration_forward.mat"); 
		return false; 
	}
	//cmask selects the centroids actually active in the microscope.
	matvar_t *cmask = NULL; 
	cmask = Mat_VarRead(matfp,"cmask"); 
	if(cmask != NULL){
		if(cmask->rank == 2 && cmask->dims[0] == 3000 && cmask->dims[1] == 1){
			if(cmask->data_type == MAT_T_UINT8){
				g_activeCentroids = 0; 
				unsigned char* dat = (unsigned char*)(cmask->data); 
				for(int i=0; i<3000; i++){
					g_cmask[i] = dat[i]; 
					if(dat[i] > 0)
						g_activeCentroids++; 
				}
				printf("found %d active centroids\n", g_activeCentroids); 
			} else {
				printf("variable cmask in calibration_forward.mat should be type unsigned char\n"); 
			}
		} else {
			printf("variable cmask in calibration_forward.mat should be 3000 x 1 \n"); 
		}
		Mat_VarFree(cmask); 
	} else {
		printf("could not load variable cmask from calibration_forward.mat\n"); 
	}
	res = read_matfile_variable(matfp, "data/calibration_forward.mat", "Cforward", g_activeCentroids*2 + 1, 97, &g_cforward); 
	Mat_Close(matfp);
	if(!res) return res; 
	
	for(int j=0; j<geneopt_num; j++){
		matfp = Mat_Open(geneopt_fnames[j], MAT_ACC_RDONLY); 
		if ( NULL == matfp ) { 
			fprintf(stderr,"Error opening %s", geneopt_fnames[j]); 
			return false; 
		}
		res = read_matfile_variable(matfp, geneopt_names[j], "genecalib", g_activeCentroids, 2, &(g_genecalib[j])); 
		if(!res){ return res; }
		Mat_Close(matfp);
	}
	// allocate a matrix for actively updated centroid info
	if(g_dat) gsl_matrix_free(g_dat); 
	g_dat = gsl_matrix_alloc(g_activeCentroids, 2); 
	
	matfp = Mat_Open("data/calibration_zernike.mat",MAT_ACC_RDONLY); 
	if ( NULL == matfp ) { 
		fprintf(stderr,"Error opening calibration_zernike.mat"); 
		return false; 
	}
	res = read_matfile_variable(matfp, "data/calibration_zernike.mat", "Z", g_activeCentroids, g_nZernike, &g_Z); 
	if(!res){ Mat_Close(matfp); return res; }
	res = read_matfile_variable(matfp, "data/calibration_zernike.mat", "dZx", g_activeCentroids, g_nZernike, &g_dZx); 
	if(!res){ Mat_Close(matfp); return res; }
	res = read_matfile_variable(matfp, "data/calibration_zernike.mat", "dZy", g_activeCentroids, g_nZernike, &g_dZy); 
	if(!res){ Mat_Close(matfp); return res; }
	Mat_Close(matfp);
	if(g_zcoef) gsl_matrix_free(g_zcoef); 
	g_zcoef = gsl_matrix_alloc(g_nZernike, 1); 
	for(int i=0; i<g_nZernike; i++){
		gsl_matrix_set(g_zcoef, i, 0, 0.0); 
	}
	
	printf("All DM control variables loaded properly.\n"); 
	return true; 
}

void dm_control_run(float* zernike, int geneopt_active, float* command)
// zernike is eg 36x1; command is eg 97x1.
{
	int k = 0; 
	for(int i=0; i<3000 && k<g_activeCentroids; i++){
		if(g_cmask[i]){
			gsl_matrix_set(g_dat, k, 0, g_centroids[i][0]); 
			gsl_matrix_set(g_dat, k, 1, g_centroids[i][1]); 
			k++; 
		}
	}
	gsl_matrix_sub(g_dat, g_genecalib[geneopt_active]); // this is dx and dy, flat-compensated. output is on g_dat.
	
	for(int i=0; i<g_nZernike; i++){
		gsl_matrix_set(g_zcoef, i, 0, zernike[i]); 
	}
#define CBNT CblasNoTrans
#define CBT CblasTrans
	gsl_matrix* desdx = gsl_matrix_alloc(g_activeCentroids, 1); 
	gsl_matrix* desdy = gsl_matrix_alloc(g_activeCentroids, 1); 
	gsl_blas_dgemm(CBNT, CBNT, 1.0, g_dZx, g_zcoef, 0.0, desdx); 
	gsl_blas_dgemm(CBNT, CBNT, 1.0, g_dZy, g_zcoef, 0.0, desdy); 
	
	//remove tip & tilt from measured dx and dy
	double sx = 0.0; 
	double sy = 0.0; 
	for(int i=0; i<g_activeCentroids; i++){
		sx += gsl_matrix_get(g_dat, i, 0); 
		sy += gsl_matrix_get(g_dat, i, 1); 
	}
	sx /= (double)g_activeCentroids; 
	sy /= (double)g_activeCentroids; 
	// and combine with the desired dx and dy to get 'A'
	gsl_matrix* A = gsl_matrix_alloc(g_activeCentroids*2+1, 1); 
	for(int i=0; i<g_activeCentroids; i++){
		gsl_matrix_set(A, i, 0,  
			gsl_matrix_get(g_dat, i, 0) - sx - gsl_matrix_get(desdx, i, 0)); 
		gsl_matrix_set(A, i + g_activeCentroids, 0,  
			gsl_matrix_get(g_dat, i, 1) - sy - gsl_matrix_get(desdy, i, 0)); 
	}
	gsl_matrix_set(A, g_activeCentroids*2, 0, 1.0); 
	
	gsl_blas_dgemm(CBT, CBNT, -0.5, g_cforward, A, 0.995, g_dmcommand);
	
	for(int i=0; i<97; i++){
		double x = gsl_matrix_get(g_dmcommand, i, 0); 
		if(x > 0.15) x = 0.15; 
		if(x<-0.15) x = -0.15; 
		command[i] = x; 
		//piston, tip and tilt will be removed on the other thread. 
	}
	gsl_matrix_free(desdx); 
	gsl_matrix_free(desdy); 
}

//init gaussian random number generator
std::default_random_engine generator;
std::normal_distribution<double> distribution(0.0,1.0);
void dm_rand_stim(float* command)
{
	//command = Z * randn(36, 1) + randn(97, 1) .* 0.022; 
	if(g_dmcommand == NULL || g_dmZ == NULL){
		printf("g_dmcommand variable not initialized .. bad file?\n"); 
		return; 
	}
	for(int i=0; i<97; i++){
		double r = distribution(generator); 
		gsl_matrix_set(g_dmcommand, i, 0, r); 
	}
	gsl_matrix* rz = gsl_matrix_alloc(g_nZernike, 1); 
	for(int i=0; i<g_nZernike; i++){
		double r = distribution(generator); 
		gsl_matrix_set(rz, i, 0, r); 
	}
	gsl_blas_dgemm(CBNT, CBNT, 1.0, g_dmZ, rz, 0.022, g_dmcommand); 
	
	for(int i=0; i<97; i++){
		double x = gsl_matrix_get(g_dmcommand, i, 0); 
		if(x > 0.15) x = 0.15; 
		if(x<-0.15) x = -0.15; 
		command[i] = x; 
		//piston, tip and tilt will be removed on the other thread. 
	}
	gsl_matrix_free(rz); 
}

void dm_control_cleanup()
{
	if(g_cforward) gsl_matrix_free(g_cforward); 
	for(int j=0; j<geneopt_num; j++){
		if(g_genecalib[j]) gsl_matrix_free(g_genecalib[j]);
	}
	if(g_dat) gsl_matrix_free(g_dat); 
	if(g_zcoef) gsl_matrix_free(g_zcoef); 
	if(g_Z) gsl_matrix_free(g_Z); 
	if(g_dZx) gsl_matrix_free(g_dZx);
	if(g_dZy) gsl_matrix_free(g_dZy);
	if(g_dmcommand) gsl_matrix_free(g_dmcommand);
	if(g_dmZ) gsl_matrix_free(g_dmZ);
}

int dm_control_geneopt_num(){
	return geneopt_num; 
}
const char* dm_control_geneopt_name(int indx){
	return geneopt_names[indx]; 
}
