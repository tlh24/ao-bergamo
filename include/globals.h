#ifndef __GLOBALS_H__
#define __GLOBALS_H__
// gl-main.cpp

extern bool g_die; 

extern unsigned char* g_data[4]; 
extern unsigned int g_w; 
extern unsigned int g_h; 
extern unsigned int g_s[4]; // in bytes. 
extern unsigned char g_copy[4]; 

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
extern Semaphore* dm_semaphore; 

void* dmControl_thread(void* writeData); 

#define MAX_LENSLETS 3000
#define NTHREADS 8
extern float g_centroids[MAX_LENSLETS][2];
extern float g_centroidsCalib[MAX_LENSLETS][2];
extern int g_lensletStarts[MAX_LENSLETS][2];
extern int g_nCentroids;
extern int g_exposure; 
extern bool g_set_exposure; 
extern bool g_calibrated; 
extern bool g_reset_data; //clear memory
extern bool g_write_data; //write to disc
extern bool g_record_data; 
extern int g_nFrames;
extern Gtk_UpdateLabel g_centroidCalc_label;
extern Gtk_UpdateLabel g_framerate_label;
extern Gtk_UpdateLabel g_dataSize_label;


void* video_thread(void*);

#endif
