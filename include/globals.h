#ifndef __GLOBALS_H__
#define __GLOBALS_H__
// gl-main.cpp

extern bool g_die; 

extern unsigned char* g_data[4]; 
extern unsigned int g_w; 
extern unsigned int g_h; 
extern unsigned int g_s[4]; // in bytes. 
extern unsigned char g_copy[4]; 

void printf_log(const char* fmt, ...);

class Gtk_UpdateLabel{
	// for displaying online metrics etc. 
private:
	std::string	m_label; 
	float		m_val; 
	GtkWidget* m_widget; 
public:	
	Gtk_UpdateLabel(){
		m_label = "none yet"; 
		m_widget = NULL; 
	}
	~Gtk_UpdateLabel(){
		if(m_widget) delete m_widget; 
	}
	void set(float v){ 
		m_val = v; 
	}
	void make(const char* lbl, GtkWidget* box){
		//make the gtk widget & add to a layout box. 
		m_label = lbl; 
		m_widget = gtk_label_new(m_label.c_str()); 
		gtk_box_pack_start (GTK_BOX (box), m_widget, FALSE, FALSE, 2);
	}
	void redraw(){
		std::string lbl = m_label + std::to_string(m_val); 
		gtk_label_set_text(GTK_LABEL(m_widget), lbl.c_str()); 
		gtk_widget_queue_draw(m_widget); 
	}
};
// in gl-main.cpp
void Gtk_CheckboxLabel_clicked(GtkWidget *widget, gpointer gdata); 

class Gtk_CheckboxLabel{
private:
	bool	m_val; 
	GtkWidget* m_widget; 
public:
	Gtk_CheckboxLabel(bool b){
		m_widget = NULL; 
		m_val = b; 
	}
	~Gtk_CheckboxLabel(){
		if(m_widget) delete m_widget; 
	}
	void set(bool b){
		if(b != m_val){
			m_val = b; 
			if(m_widget)
				gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(m_widget), b); 
		}
	}
	bool get(){ return m_val; }
	void clicked() { 
		m_val = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(m_widget));
	}
	void make(const char* lbl, GtkWidget* box){
		m_widget = gtk_check_button_new_with_label(lbl); 
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(m_widget), m_val); 
		g_signal_connect(m_widget,"clicked", 
								G_CALLBACK(Gtk_CheckboxLabel_clicked), this);
		gtk_box_pack_start (GTK_BOX (box), m_widget, FALSE, FALSE, 2);
	}
};

class Semaphore{
private:
    std::mutex mutex_;
    std::condition_variable condition_;
    unsigned long count_ = 0; // Initialized as locked.

public:
    void notify() {
        std::lock_guard<decltype(mutex_)> lock(mutex_);
        ++count_;
        condition_.notify_one();
    }

    void wait() {
        std::unique_lock<decltype(mutex_)> lock(mutex_);
        while(!count_) // Handle spurious wake-ups.
            condition_.wait(lock);
        --count_;
    }

    bool try_wait() {
        std::lock_guard<decltype(mutex_)> lock(mutex_);
        if(count_) {
            --count_;
            return true;
        }
        return false;
    }
};

#define MAX_LENSLETS 3000
#define NTHREADS 8

extern Semaphore* dm_semaphore; 
extern Semaphore g_calc_centroid_semaphore[NTHREADS]; 
extern Semaphore g_done_centroid_semaphore[NTHREADS]; 

void* dmControl_thread(void* writeData); 

extern float g_centroids[MAX_LENSLETS][2];
extern float g_centroidsCalib[MAX_LENSLETS][2];
extern int g_lensletStarts[MAX_LENSLETS][2];
extern int g_nCentroids;
extern int g_exposure; 
extern int g_actuator; 
extern bool g_set_exposure; 
extern bool g_calibrated; 
extern bool g_reset_data; //clear memory
extern bool g_write_data; //write to disc
extern Gtk_CheckboxLabel g_test_dm; 
extern Gtk_CheckboxLabel g_control_dm;
extern Gtk_CheckboxLabel g_record_data; 
extern int g_nFrames;
extern Gtk_UpdateLabel g_centroidCalc_label;
extern Gtk_UpdateLabel g_framerate_label;
extern Gtk_UpdateLabel g_dataSize_label;


void* video_thread(void*);

bool dm_control_init(); 
void dm_control_run(float* zernike, float* ctrl);
void dm_rand_stim(float* ctrl); 
void dm_control_cleanup(); 

#endif
