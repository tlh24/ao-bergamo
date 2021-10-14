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

#include <cnpy.h>
#include "gettime.h"
#include "sock.h"

const int64_t ncentroids = 1075;

using namespace torch;


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
					u[i][j] = loaded_data[i*arr.shape[1] + j]; 
				}
			}
		}
		if(arr.shape.size() == 1){
			torch::Tensor u = val.value(); 
			for(int i = 0; i < arr.shape[0]; i++){
				u[i] = loaded_data[i]; 
			}
		}
	}
}
#define CMD_SIZ 5120
static float g_cmdt[CMD_SIZ]; 
int g_controlSock = 0; //RX from matlab / ScanImage.

void my_handler(int s){
	printf("Caught signal %d\n",s);
	if(g_controlSock)
		close_socket(g_controlSock); 
	exit(1); 
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
	if (torch::cuda::is_available()) {
		std::cout << "CUDA is available! Controller on GPU." << std::endl;
		device = torch::Device(torch::kCUDA);
	}
	// note: both GPU and CPU have about the same latency on the Lenovo Carbon X1-v3 - Geforce 1650

	dmControlNet controller(ncentroids); 
	LoadStateDict(controller, "/home/tlh24/ao-bergamo/ml/dmcontrolnet"); 
	controller->to(device); 
	for(int i=0; i<10; i++){
		long double sta = gettime(); 
		torch::Tensor random_wavefront = torch::randn({97});
		random_wavefront = random_wavefront.to(device); 
		torch::Tensor dmctrl = controller->forward(random_wavefront); 
		std::cout << "computation time " << (gettime() - sta)*1000.0 << " ms" << std::endl; 
		//without transferring random_wavefront to the gpu, computation is ~130us (!!)
	}
	// std::cout << dmctrl << std::endl; 
	// this should be ~= BestDMcommand.
	// and indeed it looks ok (not perfect?)

	g_controlSock = setup_socket(13131); 
	std::cout << "UDP socket listening on port 13131" << std::endl; 
	while(1){
		bzero(g_cmdt, sizeof(g_cmdt)); 
		int n = recvfrom(g_controlSock, (char*)g_cmdt, CMD_SIZ,0,0,0); 
		if( n > 0 ){
			//std::cout << "got " << n << " data bytes" << std::endl; 
			if( n == 98 * sizeof(float) ){
				float check = g_cmdt[0]; 
				if(check > 3.141 && check < 3.1416){
					long double sta = gettime(); 
					torch::Tensor wf = torch::from_blob(&(g_cmdt[1]), {97}); 
					wf = wf.to(device); 
					torch::Tensor dmctrl = controller->forward(wf);
					dmctrl = dmctrl.to(torch::kCPU); 
					// need to send this to the deformable mirror..
					std::cout << "computation time " << (gettime() - sta)*1000.0 << " ms" << std::endl; 
					//std::cout << dmctrl << std::endl;
				}
			}
		}
	}
}
