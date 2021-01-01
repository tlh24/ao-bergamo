// Alpao SDK Header: All types and class are in the ACS namespace
#include "asdkDM.h"
using namespace acs;

// System Headers
#include <iostream>
#include <unistd.h>
using namespace std;

#include <mutex>
#include <condition_variable>

#include <gtk/gtk.h>
#include <gdk/gdk.h>

#include "globals.h"

Semaphore* dm_semaphore; 

// Example using C++ API
void* dmControl_thread(void* writeData)
{
	float* fdata = (float*)writeData; 
	String serial;
	serial = "BAX390"; 

	// Load configuration file
	DM dm( serial.c_str() );

	// Get the number of actuators
	UInt nbAct = (UInt) dm.Get( "NbOfActuator" );

	// Check errors
	if ( !dm.Check() ){
		return NULL;
	}
	cout << "Number of actuators: " << nbAct << endl;
	
	dm_semaphore = new Semaphore(); 

	// Initialize data
	Scalar *data = new Scalar[nbAct];
	for(int i=0; i<nbAct; i++){
		data[i] = fdata[i]; 
	}

	while(!g_die){
		dm_semaphore->wait(); 
		cout << "Sending data to mirrors" << endl;
		for(int i=0; i<nbAct; i++)
			data[i] = fdata[i]; 
		// Send value to the DM
		dm.Send( data );
		cout << "Done." << endl;
		DM::PrintLastError();
	}
	// Release memory
	delete [] data;

	return NULL;
}
