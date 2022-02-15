#include <torch/torch.h>
#include <torch/csrc/jit/serialization/import.h>

#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_fit.h>

#include <cnpy.h>
#include "gettime.h"
#include "sock.h"

#define ALPAO 1
#ifdef ALPAO
// Alpao SDK Header: All types and class are in the ACS namespace
#include "asdkDM.h"
#endif

const int64_t ncentroids = 1075;

using namespace torch;
using namespace torch::indexing;


struct dmControlNetImpl : nn::Module {
	dmControlNetImpl(int ncentroids) : 
		linearVS(nn::LinearOptions(97, ncentroids)), 
		linear1(nn::LinearOptions(ncentroids, 1024)), 
		relu1(nn::LeakyReLUOptions().negative_slope(0.2)), 
		linear2(nn::LinearOptions(1024, 512)), 
		relu2(nn::LeakyReLUOptions().negative_slope(0.2)), 
		linear3(nn::LinearOptions(512, 256)), 
		relu3(nn::LeakyReLUOptions().negative_slope(0.2)), 
		linear4(nn::LinearOptions(256, 97))
	{
		// register_module() is needed if we want to use the parameters() method later on
		// first layer uses 'bias' to set the flat wavefront
		register_module("VS", linearVS); 
		register_module("0", linear1); 
		register_module("2", linear2); 
		register_module("4", linear3); 
		register_module("6", linear4); 
		register_module("relu1", relu1); 
		register_module("relu2", relu2); 
		register_module("relu3", relu3); 
	}

	torch::Tensor forward(torch::Tensor x) {
		x = linearVS(x);
		x = relu1(linear1(x));
		x = relu2(linear2(x));
		x = relu3(linear3(x));
		x = linear4(x);
		return x * 0.15; //NB! training scaled to +-1.0
	}
	torch::Tensor vs_slice(int i){
		torch::Tensor vs = linearVS->weight; 
		return vs.index({Slice(), i}); 
	}

	nn::Linear linearVS, linear1, linear2, linear3, linear4;
	nn::LeakyReLU relu1, relu2, relu3; 
};
TORCH_MODULE(dmControlNet);


void LoadStateDict(dmControlNet& module,
							const std::string& dir_name) {
	torch::NoGradGuard no_grad;
	auto params = module->named_parameters(true /*recurse*/);
	for (auto& val : params) {
		std::string file_name = dir_name + "/" + val.key() + ".npy"; 
		cnpy::NpyArray arr = cnpy::npy_load(file_name); 
		float* loaded_data = arr.data<float>(); 
		std::cout << "loaded " << file_name << " shape " << arr.shape << std::endl; 
		if(arr.shape.size() == 2){
			torch::Tensor u = val.value(); 
			for(int i = 0; i < arr.shape[0]; i++){
				for(int j = 0; j < arr.shape[1]; j++){
					u.index_put_({i,j}, loaded_data[i*arr.shape[1] + j]); 
				}
			}
		}
		if(arr.shape.size() == 1){
			torch::Tensor u = val.value(); 
			for(int i = 0; i < arr.shape[0]; i++){
				u.index_put_({i}, loaded_data[i]); 
			}
		}
	}
}
#define CMD_SIZ 5120
static float g_cmdt[CMD_SIZ]; 
int g_controlSock = 0; //RX from matlab / ScanImage.
bool g_die = false; 

void my_handler(int s){
	printf("Caught signal %d\n",s);
	g_die = true;  
}

void dm_actuator_to_xy(int act, double* x, double* y){
	int starts[] = {0,5,12,21,32,43,54,65,76,85,92,  97};
	int startsY[] = {-2,-3,-4,-5,-5,-5,-5,-5,-4,-3,-2}; 
	int startsX[] = {5,4,3,2,1,0,-1,-2,-3,-4,-5}; 
	if(act < 0 || act > 96) return; 
	int i = 0; 
	while(i<11 && starts[i+1] <= act){
		i++; 
	}
	*x = (double)( startsX[i] ); 
	*y = (double)( startsY[i] + act - starts[i] ); 
}

void dm_remove_ptt(double* dm_data){
	// remove piston, tilt, and tip from the DM control signal.  
	float piston = 0.f; 
	for(int i=0; i<97; i++){
		piston += dm_data[i]; 
	}
	piston /= 97.f; 
	for(int i=0; i<97; i++){
		dm_data[i] -= piston; 
	}
	double x[97]; //gsl routines are double prec
	double y[97]; 
	double dm[97]; 
	for(int i=0; i<97; i++){
		dm_actuator_to_xy(i, &(x[i]), &(y[i]));
		dm[i] = dm_data[i]; 
	}
	// use X and Y to predict dm via linear regression
	double cx, cov11, sumsq; 
	gsl_fit_mul(x, 1, dm, 1, 97, &cx, &cov11, &sumsq); 
	for(int i=0; i<97; i++){
		dm[i] = dm[i] - x[i]*cx; 
	}
	double cy; 
	gsl_fit_mul(y, 1, dm, 1, 97, &cy, &cov11, &sumsq); 
	for(int i=0; i<97; i++){
		dm_data[i] = dm[i] - y[i]*cy; 
	}
}
	
int main(int argc, const char* argv[]) {
	g_startTime = gettime(); 
	torch::manual_seed(1);
	
	// capture control-C to terminate program / close socket. 
	struct sigaction sigIntHandler;
	sigIntHandler.sa_handler = my_handler;
	sigemptyset(&sigIntHandler.sa_mask);
	sigIntHandler.sa_flags = 0;
	sigaction(SIGINT, &sigIntHandler, NULL);

	// Create the device we pass around based on whether CUDA is available.
	torch::Device device(torch::kCPU);
// 	if (torch::cuda::is_available() && 0) {
// 		std::cout << "CUDA is available! Controller on GPU." << std::endl;
// 		device = torch::Device(torch::kCUDA);
// 	}
	// note: both GPU and CPU have about the same latency on the Lenovo Carbon X1-v3 - Geforce 1650

	dmControlNet controller(ncentroids); 
	LoadStateDict(controller, "/home/tlh24/ao-bergamo/ml/dmcontrolnet"); 
	controller->to(device); 
	torch::Tensor dmctrl;
	for(int i=0; i<10; i++){
		long double sta = gettime(); 
		torch::Tensor random_wavefront = torch::zeros({97});
		random_wavefront = random_wavefront.to(device); 
		dmctrl = controller->forward(random_wavefront); 
		std::cout << "computation time " << (gettime() - sta)*1000.0 << " ms" << std::endl; 
		//without transferring random_wavefront to the gpu, computation is ~130us (!!)
	}
	std::cout << dmctrl << std::endl; 
	// this should be ~= BestDMcommand.
	// and indeed it looks ok (not perfect?)
	
#ifdef ALPAO
	std::string serial = "AlpaoLib/BAX390"; 
	// Load configuration file
	acs::DM dm( serial.c_str() );

	// Get the number of actuators
	acs::UInt nbAct = (acs::UInt) dm.Get( "NbOfActuator" );

	// Check errors
	if ( !dm.Check() ){
		std::cout << "deformable mirror is not responding" << std::endl; 
		return 1;
	}
	std::cout << "Number of actuators: " << nbAct << std::endl;

	// Initialize data
	acs::Scalar *data = new acs::Scalar[nbAct];
	for(int i=0; i<nbAct; i++){
		data[i] = 0.0; 
	}
	
// 	std::cout << controller->vs_slice(0) << std::endl; 

	g_controlSock = setup_socket(13131); 
	std::cout << "UDP socket listening on port 13131" << std::endl; 
	while(!g_die){
		bzero(g_cmdt, sizeof(g_cmdt)); 
		int n = recvfrom(g_controlSock, (char*)g_cmdt, CMD_SIZ,0,0,0); 
		if( n > 0 ){
			//std::cout << "got " << n << " data bytes" << std::endl; 
			if( n == 98 * sizeof(float) ){
				float check = g_cmdt[0]; 
				if(check > 3.141 && check < 3.1416){
					long double sta = gettime(); 
					torch::Tensor wf = torch::from_blob(&(g_cmdt[1]), {97}); 
// 					torch::Tensor wf = torch::zeros({97}); 
// 					wf[3] = 14.0*sin((double)gettime()); 
					wf = wf.to(device); 
					torch::Tensor dmctrl = controller->forward(wf);
					dmctrl = dmctrl.to(torch::kCPU); 
					long double comp = gettime(); 
					float* p = dmctrl.data_ptr<float>();
					for(int i=0; i<nbAct; i++){
						data[i] = p[i]; 
						//clamp values to keep from damaging the mirror.
						if(data[i] > 0.15) data[i] = 0.15; 
						if(data[i] < -0.15) data[i] = -0.15; 
					}
					dm_remove_ptt( data ); 
					dm.Send( data );
// 				std::cout << dmctrl << std::endl; 
					std::cout << "computation " << (comp - sta)*1000.0 << " ms "; 
					std::cout << "total " << (gettime() - sta)*1000.0 << " ms\n"; 
				}
				if(check > 2.718 && check < 2.719){
					long double sta = gettime(); 
					for(int i=0; i<97; i++){
						float a = g_cmdt[i+1] * 0.15; 
						if(a > 0.15) a = 0.15; 
						if(a <-0.15) a = -0.15; 
						data[i] = a; 
					}
					dm_remove_ptt( data ); 
					dm.Send( data );
					std::cout << "total " << (gettime() - sta)*1000.0 << " ms\n"; 
				}
			}
		}
	}
	delete [] data;
	acs::DM::PrintLastError();
	close_socket(g_controlSock); 
#endif
}
