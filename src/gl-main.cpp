#include	<stdio.h>
#include	<stdlib.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <string>
#include <mutex>
#include <condition_variable>
#include <gsl/gsl_linalg.h>

//opengl headers.
#include <gtk/gtk.h>
#include <gdk/gdk.h>
// #include <gdk/gdkx.h> -- conflict with protobuf
#include <GL/glew.h>
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glx.h> 
#include <GL/glext.h>
#include "text_helper.h"

#include "gettime.h"
#include "globals.h"

// Globals
bool g_die = false; 
bool g_render_targets = true;

GtkWidget	*glarea = NULL;
GtkWidget	*g_logarea = NULL; 
Gtk_UpdateLabel g_centroidCalc_label;
Gtk_UpdateLabel g_framerate_label;
Gtk_UpdateLabel g_dataSize_label;

GLuint	gl_tex[4]; 
GLuint   gl_pbo[4][2]; 
int gl_pbo_indx[4] = {0,0,0,0}; 

unsigned char* g_data[4]; 
unsigned int g_w; 
unsigned int g_h; 
unsigned int g_s[4]; // GL size, in bytes. 
unsigned char g_copy[4]; 
int g_viewport_width = 2048; 
int g_viewport_height = 2048; 
unsigned int	g_frame; // OpenGL frame number!

//Text and font stuff
struct point {
    GLfloat x;
    GLfloat y;
    GLfloat s;
    GLfloat t;
};

GLuint program;
GLint attribute_coord;
GLint uniform_tex;
GLint uniform_color;

#define clamp(a,b,c)( a < b ? b : (a > c ? c : a) )

// model-view-projection matrix
float mvp[16];

// GL objects 
guint g_vbo[6] = {0,0,0,0,0,0};
guint g_vao = 0; 
guint g_program[3] = {0,0,0};
guint g_mvp_location[3] = {0,0,0};
guint g_flip_location = 0;
guint g_texsamp_location; 
guint g_position_location[3];
guint g_tex_location;
guint g_color_location;
float	g_mvp[16]; 

/* the vertex data is constant (a square) */
static const float vertex_data[] = {
   -1.f,  1.f, 
	-1.f, -1.f, 
	 1.f,  1.f, 
	 1.f, -1.f, 
	-1.f, -1.f, 
	 1.f,  1.f
};

static const float uv_data[] = {
	0.f, 0.f, 
	0.f, 1.f, 
	1.f, 0.f, 
	1.f, 1.f, 
	0.f, 1.f, 
	1.f, 0.f
};

static const float indicator_data[] = {
	0.f, 1.f, //ccw draw
	0.07f, 1.f, 
	0.f, 0.f, 
	1.f, 0.f, 
	1.f, -0.07f, 
	0.f, 0.f, 
	0.f, -1.f, 
	-0.07f, -1.f, 
	0.f, 0.f, 
	-1.f, 0.f, 
	-1.f, 0.07f, 
	0.f, 0.f
};

GtkWidget *g_exposureSpin[1] = {0};
bool g_display_camimage = true; 
bool g_display_centroids = true; 

int init_resources() {
	/* Initialize the FreeType2 library */
	printf("init text stuff \n");
	init_text_helper(g_vao, g_vbo[4]);
	return 1;
}

void printf_log(const char* fmt, ...){
	if(g_logarea){
		GtkTextBuffer* txtBuf = 
			gtk_text_view_get_buffer( GTK_TEXT_VIEW(g_logarea)); 
		GtkTextIter iter; 
		gtk_text_buffer_get_end_iter(txtBuf, &iter);
		char str[256]; 
		va_list ar; 
		va_start(ar, fmt);
		vsnprintf(str, 256, fmt, ar);
		va_end(ar);
		gtk_text_buffer_insert( txtBuf, &iter, str, -1); 
	}
}

gint destroyGUI(GtkWidget *widget, gpointer gdata){
	g_die = true; 
	sleep(1); //wait for video threads to close. 
	gtk_main_quit();
	return (FALSE);
}

char* read_txt_file(const char* fname){
	char *source = NULL;
	FILE *fp = fopen(fname, "r");
	if (fp != NULL) {
		/* Go to the end of the file. */
		if (fseek(fp, 0L, SEEK_END) == 0) {
			/* Get the size of the file. */
			long bufsize = ftell(fp);
			if (bufsize == -1) { /* Error */ }
			/* Allocate our buffer to that size. */
			source = (char*)malloc(sizeof(char) * (bufsize + 1));
			/* Go back to the start of the file. */
			if (fseek(fp, 0L, SEEK_SET) != 0) { /* Error */ }
			/* Read the entire file into memory. */
			size_t newLen = fread(source, sizeof(char), bufsize, fp);
			if (newLen == 0) {
					fputs("Error reading file", stderr);
			} else {
					source[newLen++] = '\0'; /* Just to be safe. */
			}
		}
		fclose(fp);
		return source; 
	}else
		return 0; 
}

static guint create_shader (int          shader_type,
               const char  *source,
               guint       *shader_out)
{
	guint shader = glCreateShader (shader_type);
	glShaderSource (shader, 1, &source, NULL);
	glCompileShader (shader);

	int status;
	glGetShaderiv (shader, GL_COMPILE_STATUS, &status);
	if (status == GL_FALSE){
		int log_len;
		glGetShaderiv (shader, GL_INFO_LOG_LENGTH, &log_len);

		char *buffer = (char*)g_malloc (log_len + 1);
		glGetShaderInfoLog (shader, log_len, NULL, buffer);

		printf("Compilation failure in %s shader: %s",
			shader_type == GL_VERTEX_SHADER ? "vertex" : "fragment",
			buffer);

		g_free (buffer);

		glDeleteShader (shader);
		shader = 0;
	}

	if (shader_out != NULL)
		*shader_out = shader;

	return shader != 0;
}

void setCentroidData(){
	glBindBuffer(GL_ARRAY_BUFFER, g_vbo[3]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(g_centroids), g_centroids, GL_DYNAMIC_DRAW); 
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, (void*)0); //XY 
	//this call is critical!!
}

static gboolean init_shaders () {
	char* src;
	guint vertex = 0, fragment = 0;
	int status = 0;

	for(int i=0; i<3; i++){
		/* load the vertex shader */
		if(i == 0)
			src = read_txt_file("gl/vertex.glsl");
		else if(i == 1)
			src = read_txt_file("gl/vertex_flat.glsl");
		else if(i == 2)
			src = read_txt_file("gl/text.v.glsl");

		if(src){
			create_shader (GL_VERTEX_SHADER, src, &vertex);
			free(src); 
		}

		/* load the fragment shader */
		if(i == 0)
			src = read_txt_file("gl/fragment.glsl"); 
		else if(i == 1)
			src = read_txt_file("gl/fragment_flat.glsl"); 
		else if(i == 2)
			src = read_txt_file("gl/text.f.glsl");
		if(src){
			create_shader (GL_FRAGMENT_SHADER, src, &fragment);
			free(src); 
		}
		std::cout<<vertex<<std::endl;
		std::cout<<fragment<<std::endl;

		/* link the vertex and fragment shaders together */
		if(vertex && fragment){
			g_program[i] = glCreateProgram ();
			glAttachShader (g_program[i], vertex);
			glAttachShader (g_program[i], fragment);
			glLinkProgram (g_program[i]);

			glGetProgramiv (g_program[i], GL_LINK_STATUS, &status);
			if (status == GL_FALSE){
				int log_len = 0;
				glGetProgramiv (g_program[i], GL_INFO_LOG_LENGTH, &log_len);

				char *buffer = (char*)g_malloc(log_len + 1);
				glGetProgramInfoLog (g_program[i], log_len, NULL, buffer);

				printf_log("Linking failure in program: %s", buffer);

				g_free (buffer);
				glDeleteProgram (g_program[i]);
				g_program[i] = 0;
			}

			/* get the location of the "mvp" uniform */
			g_mvp_location[i] = glGetUniformLocation (g_program[i], "mvp");
			printf_log("g_mvp_location %d\n", g_mvp_location[i]); 
			
			if(i==0){
				g_texsamp_location = glGetUniformLocation(g_program[i], "myTexture");
			
				/* get the location of the "position" and "color" attributes */
				g_position_location[i] = glGetAttribLocation(g_program[i], "position");
				g_tex_location = glGetAttribLocation(g_program[i], "vertTexCoord");
				g_flip_location = glGetUniformLocation(g_program[i], "flip");
				printf_log("position %d, vertTexCoord %d flip %d\n", g_position_location[1], g_tex_location, g_flip_location); 
			}
			if(i==1){
				g_color_location = glGetUniformLocation(g_program[i], "myColor");
				g_position_location[i] = glGetAttribLocation(g_program[i], "position");
			}
			if(i==2){
			}
			/* the individual shaders can be detached and destroyed */
			glDetachShader (g_program[i], vertex);
			glDetachShader (g_program[i], fragment);
		}
	}
	glUseProgram(g_program[1]);
	
	glGenVertexArrays(1, &g_vao);
	glBindVertexArray(g_vao);
	
	//init the vertex buffer objects. 
	glGenBuffers(6, &(g_vbo[0])); 
	glBindBuffer(GL_ARRAY_BUFFER, g_vbo[0]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertex_data), vertex_data, GL_STATIC_DRAW);
	
	glBindBuffer(GL_ARRAY_BUFFER, g_vbo[1]); 
	glBufferData(GL_ARRAY_BUFFER, sizeof(uv_data), uv_data, GL_STATIC_DRAW);
	
	float ind_data[24 * 5]; 
	float offsets[] = {
		0.0, 0.0, 
		0.025, 0.0, 
		-0.025, 0.0,
		0.0, 0.025, 
		0.0, -0.025 }; // ghetto anti-aliasing.
	for(int i=0; i<5; i++){
		float x = offsets[i*2 + 0]; 
		float y = offsets[i*2 + 1]; 
		for(int j=0; j<12; j++){
			ind_data[i*24 + j*2 + 0] = indicator_data[j*2 + 0] + x; 
			ind_data[i*24 + j*2 + 1] = indicator_data[j*2 + 1] + y; 
		}
	}
	glBindBuffer(GL_ARRAY_BUFFER, g_vbo[2]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(ind_data), ind_data, GL_STATIC_DRAW);
  
	float circle_data[12*3*2]; 
	//12 triangles, radial. 12 * 3 * 2 = 72. 
	float t = 0.0; 
	float dt = 2.0 * 3.1415926 / 12.0; 
	for(int i=0; i<12; i++){
		circle_data[i*6+0] = cos(t); 
		circle_data[i*6+1] = sin(t); 
		circle_data[i*6+2] = cos(t+dt); 
		circle_data[i*6+3] = sin(t+dt); 
		circle_data[i*6+4] = 0.0; 
		circle_data[i*6+5] = 0.0; 
		t += dt; 
	}
	glBindBuffer(GL_ARRAY_BUFFER, g_vbo[5]); //text is #4.
	glBufferData(GL_ARRAY_BUFFER, sizeof(circle_data), circle_data, GL_STATIC_DRAW);

	glBindVertexArray(0); //release vao
	
	return g_program[0] != 0;
}

//matrix math -- should probably be a library, but easy to write!
//matrices are C-style, row-major. 
void mat_identity(float* m){
	for(int i=0; i<16; i++)
		m[i] = 0.f; 
	for(int i=0; i<4; i++)
		m[i*5] = 1.f; 
}
void mat_mul(float* a, float* b, float *d){
	//d = a * b
	for(int i=0; i<16; i++)
		d[i] = 0.f; 
	for(int r=0; r<4; r++){
		for(int c=0; c<4; c++){
			for(int j=0; j<4; j++)
				d[r*4+c] += a[r*4+j] * b[j*4+c]; 
		}
	}
}
void mat_vec(float* a, float* b, float *d){
	// d = a * b, where b is a column vector. 
	for(int i=0; i<4; i++)
		d[i] = 0.f; 
	for(int r=0; r<4; r++){
		for(int j=0; j<4; j++)
			d[r] += a[j*4+r] * b[j]; //note the transpose.
	}
}
void mat_inverse(float* a, float* b){
	double ad[16]; 
	double bd[16]; 
	int i, s; 
	
	for(i=0; i<16; i++)
		ad[i] = a[i]; 
	
	gsl_matrix_view m = gsl_matrix_view_array(ad, 4, 4);
	gsl_matrix_view inv = gsl_matrix_view_array(bd, 4, 4);
	gsl_permutation * p = gsl_permutation_alloc(4);

	gsl_linalg_LU_decomp (&m.matrix, p, &s);    
	gsl_linalg_LU_invert (&m.matrix, p, &inv.matrix);

	gsl_permutation_free (p);
	for(i=0; i<16; i++)
		b[i] = bd[i]; 
}

void calcAR(float& xs, float& ys){
	//calculate the aspect ratio of the GL viewport.
	float ar = (float)g_viewport_width / (float)g_viewport_height; 
	xs = 1.0; ys = 1.0; 
	if(ar > 1.0){ //wider, scale x accordingly
		xs = 1.0 / ar; 
	} else { //taller, scale y. 
		ys = ar / 1.0; 
	}
}

void glErr(int ln){
	GLenum err = GL_NO_ERROR;
	while((err = glGetError()) != GL_NO_ERROR){
		printf("%d glError 0x%x\n", ln, err); 
	}
}

double last_render_time; 
static gboolean
render(GtkGLArea *area, GdkGLContext *context){
	// inside this function it's safe to use GL; the given
	// #GdkGLContext has been made current to the drawable
	// surface used by the #GtkGLArea and the viewport has
	// already been set to be the size of the allocation
	double t = gettime(); 
	GtkAllocation a; 
	gtk_widget_get_allocation((GtkWidget*)area, &a); 
	g_viewport_width = a.width; 
	g_viewport_height = a.height; 
	//printf("%d by %d\n", g_viewport_width, g_viewport_height); 
	// we can start by clearing the buffer
	glClearColor (sin(gettime()*0.8456143)*0.03f + 0.03f, 0.0, sin(gettime())*0.04f + 0.04f, 1.0);
	glErr(130879);
	glClear(GL_COLOR_BUFFER_BIT);
	glClearDepth(1);
	glDisable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	glDisable(GL_CULL_FACE); 
	glErr(130);
	//copy over the texture data. 
	for(int i=0; i<1; i++){
		if(g_display_camimage && g_copy[i] && g_data[i]){
			glBindBufferARB(GL_PIXEL_UNPACK_BUFFER_ARB, gl_pbo[i][gl_pbo_indx[i]]);
			GLubyte* ptr = (GLubyte*)glMapBufferRange(GL_PIXEL_UNPACK_BUFFER_ARB, 0, g_w*g_h, GL_MAP_WRITE_BIT|GL_MAP_INVALIDATE_BUFFER_BIT);
			if(ptr){
				memcpy((void*)ptr, g_data[i], g_s[i]); 
				glUnmapBufferARB(GL_PIXEL_UNPACK_BUFFER_ARB); // release pointer to mapping buffer
			}
			g_copy[i] = 0; 
			glBindTexture(GL_TEXTURE_2D, gl_tex[i]);
			glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, gl_pbo[i][gl_pbo_indx[i]]);
			glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, 
									g_w, g_h, GL_RED, GL_UNSIGNED_BYTE, 0);
			glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB,0);
			glBindTexture(GL_TEXTURE_2D,0);
			//gl_pbo_indx[i] = gl_pbo_indx[i]^0x1; --tried this, doesn't speed things up.
		}
	}
	glErr(845);
	
	float xs, ys; 
	calcAR(xs, ys); 
	
	float mvp[16]; 
	float flip[4]; 
	mat_identity(mvp); 
	mvp[0] = xs * 0.5; 
	mvp[5] = ys * 0.5; 
	
// 	float z[16]; 
// 	calcCamZoom(z, xs, ys); 
// 	mat_mul(mvp, z, g_mvp); //inverse used for mouse-clicks
	memcpy(g_mvp, mvp, 16*sizeof(float)); 
	
	if(g_program[0] && g_vbo[0] && g_program[1] && g_vbo[3]){
		glUseProgram (g_program[0]);
		glBindVertexArray(g_vao);
		
		glBindBuffer(GL_ARRAY_BUFFER, g_vbo[0]);
		glVertexAttribPointer(g_position_location[0], 2, GL_FLOAT, GL_FALSE, 0, (void*)0); //XY 
		glEnableVertexAttribArray(g_position_location[0]); //attribute '0' must match the layout in the shader. 
		
		glBindBuffer(GL_ARRAY_BUFFER, g_vbo[1]);
		glEnableVertexAttribArray(g_tex_location);
		glVertexAttribPointer(g_tex_location, 2, GL_FLOAT, GL_FALSE, 0, (void*)0); //UV 

		mat_identity(mvp);
		int i=0; 
		mvp[12] = 0.0; //x center
		mvp[13] = 0.0; //y center
		mvp[0] = 2.0; //x scale
		mvp[5] = 2.0; //y scale
		float md[16]; 
		mat_mul(mvp, g_mvp, md); 
		glUniformMatrix4fv(g_mvp_location[0], 1, GL_FALSE, md);
		flip[0] = 1.0; 
		flip[1] = 1.0; //flip y. 
		
		if(g_display_camimage){
			glUniform4fv(g_flip_location, 1, flip); 
			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, gl_tex[i]);
			glUniform1i(g_texsamp_location, 0);
			glDrawArrays (GL_TRIANGLES, 0, 6);
		}

		if(g_display_centroids){
			//render point sprites for the centroids. 
			setCentroidData(); 
			glUseProgram (g_program[1]);
			glUniformMatrix4fv (g_mvp_location[1], 1, GL_FALSE, md);
			glBindBuffer(GL_ARRAY_BUFFER, g_vbo[3]);
			glVertexAttribPointer(g_position_location[1], 2, GL_FLOAT, GL_FALSE, 0, (void*)0); //XY 
			glEnableVertexAttribArray(g_position_location[1]); //attribute '0' must match the layout in the shader. 
			float color[] = {1.0, 0.0, 0.0, 0.25}; 
			glUniform4fv(g_color_location, 1, color); 
			glPointSize(8.0);              //specify size of points in pixels
			glDrawArrays(GL_POINTS, 0, MAX_LENSLETS);  //draw the points
		}
	}
	glErr(602);
	glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0); // vao
	glUseProgram(0);
	glFlush ();
//   printf("render time %f copy time %f frame rate %f\n", t2 - t, t3 - t, 1.0 / (t - last_render_time)); 
  last_render_time = t; 
	return TRUE;
}

static gboolean on_configure(GtkWidget *widget, GdkEventExpose *event){
	GtkAllocation a; 
	gtk_widget_get_allocation(widget, &a); 
	g_viewport_width = a.width; 
	g_viewport_height = a.height; 
	printf("Width: %i\nHeight%i\n", a.width, a.height);
	return FALSE; 
}

static void 
on_realize (GtkGLArea *area){
	// We need to make the context current if we want to
	// call GL API
	gtk_gl_area_make_current (area);
	
	//init the GL extensions.  *Needed* for the shaders! 
	glewExperimental = true; // Needed for vao.
	GLenum err = glewInit();
	if (err != GLEW_OK)
		exit(1); // or handle the error in a nicer way
	if (!GLEW_VERSION_3_0)  // check that the machine supports the 3.0 API.
		exit(1); // or handle the error in a nicer way
	
	// If there were errors during the initialization or
	// when trying to make the context current, this
	// function will return a #GError for you to catch
	if (gtk_gl_area_get_error (area) != NULL)
		return;
		/* ... Specify clear color (blue) ... */
	glClearColor(0.0,0.0,1.0,1.0);

	glEnable(GL_TEXTURE_2D);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	//glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	glDisable(GL_CULL_FACE); 
	glEnable(GL_POINT_SPRITE);

	//glPixelStorei (GL_UNPACK_ALIGNMENT, 1); 
  
	for(int i=0; i<1; i++){
		glGenTextures(1, &(gl_tex[i]));					// Create 1 Texture
		glBindTexture(GL_TEXTURE_2D, gl_tex[i]);			// Bind The Texture
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, g_w, g_h, 0,
			GL_RED, GL_UNSIGNED_BYTE, (void*)g_data[i]);		// Build Texture Using Information In data
		glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
    //make the pbo
    glGenBuffersARB(2, gl_pbo[i]);
    for(int j=0; j<2; j++){
      glBindBufferARB(GL_PIXEL_UNPACK_BUFFER_ARB, gl_pbo[i][j]);
      glBufferDataARB(GL_PIXEL_UNPACK_BUFFER_ARB, 
                    g_w*g_h*1, 0, GL_STREAM_DRAW_ARB);
      glBindBufferARB(GL_PIXEL_UNPACK_BUFFER_ARB, 0);
    }
	}
	init_shaders ();
	init_resources();
	glErr(111);
}

static gboolean redraw(gpointer user_data){
	if(glarea){
		gtk_widget_queue_draw(glarea);
	}
	g_centroidCalc_label.redraw(); 
	g_framerate_label.redraw(); 
	g_dataSize_label.redraw(); 
	g_frame++; 
	return TRUE;
}

static void
gl_fini(){
	/* we need to ensure that the GdkGLContext is set before calling GL API */
	gtk_gl_area_make_current (GTK_GL_AREA (glarea));

	/* skip everything if we're in error state */
	if (gtk_gl_area_get_error (GTK_GL_AREA (glarea)) != NULL)
		return;

	/* destroy all the resources we created */
	if (g_vbo[0] != 0)
		glDeleteBuffers (1, g_vbo);
	for(int i=0; i<3; i++){
			if (g_program[i] != 0)
				glDeleteProgram (g_program[i]);
	}
	for(int i=0; i<4; i++){
		if(gl_tex[i]) glDeleteTextures(1, &(gl_tex[i])); 
		if(gl_pbo[i]) glDeleteBuffersARB(2, gl_pbo[i]);
		gl_tex[i] = 0; 
		if(g_data[i]) free(g_data[i]); 
		g_data[i] = 0; 
	}
}

static void updateExposureCB(GtkWidget *spinner, gpointer p){
	int k = (int)((long long)p & 0xf);
	if(k >= 0 && k < 3){
		int v = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(g_exposureSpin[k]));
		g_exposure = v; 
		g_set_exposure = true; 
	}
}

void recalibrate_clicked(GtkWidget *widget, gpointer gdata){
	g_nFrames = 0; 
	g_calibrated = false; 
	g_reset_data = true; 
}

void reset_clicked(GtkWidget *widget, gpointer gdata){
	g_reset_data = true; 
}

void write_data_clicked(GtkWidget *widget, gpointer gdata){
	g_write_data = true; 
}

void record_data_clicked(GtkWidget *widget, gpointer gdata){
	g_record_data = gtk_toggle_button_get_active(
		GTK_TOGGLE_BUTTON(widget));
}

void display_centroids_clicked(GtkWidget *widget, gpointer gdata){
	g_display_centroids = gtk_toggle_button_get_active(
		GTK_TOGGLE_BUTTON(widget));
}

void display_camimage_clicked(GtkWidget *widget, gpointer gdata){
	g_display_camimage = gtk_toggle_button_get_active(
		GTK_TOGGLE_BUTTON(widget));
}

gint glarea_motion_notify(GtkWidget *widget, GdkEventMotion *event){
	gint			x, y;
	float ex = event->x; 
	float ey = event->y; 
	GdkModifierType		state;
	GdkDevice *device = gdk_event_get_device ((GdkEvent*) event);

	if (event->is_hint) {
		gdk_window_get_device_position(event->window, device, &x, &y, &state);
	} else {
		x = (gint)event->x;
		y = (gint)event->y;
		state = (GdkModifierType)event->state;
	}

	if((state & GDK_BUTTON1_MASK) && (state & GDK_CONTROL_MASK)) {
		g_print("Mouse motion button 1 at coordinates (%d,%d)\n",x,y);
	}
	if(state & GDK_BUTTON2_MASK) {
		/* ... Zooming drag ... */
		g_print("Mouse motion button 2 at coordinates (%d,%d)\n",x,y);
	}
	if(state & GDK_BUTTON3_MASK) {
		/* ... 3rd button drag ... */
		g_print("Mouse motion button 3 at coordinates (%d,%d)\n",x,y);
	}
	return TRUE;
}

gint glarea_button_press(GtkWidget *widget, GdkEventButton *event){
	gint	return_status=TRUE;

	float cx = g_viewport_width / 5.f; 
	float cx2 = 3.0 * cx; 
	cx *= 2.f; 
	float cy = g_viewport_height / 2; 
	float ex = event->x; 
	float ey = event->y; 
	float x, y; 
	GdkModifierType state = (GdkModifierType)event->state;
	switch(event->type) {
	case GDK_BUTTON_PRESS:
		break;
	case GDK_2BUTTON_PRESS:
		break;
	case GDK_3BUTTON_PRESS:
		g_print("Mouse button %d triple-click at coordinates (%lf,%lf)\n",
				event->button,event->x,event->y);
		break;
	default:
		g_print("Unknown button press event\n");
		return_status=FALSE;
		break;
	}
	return(return_status);
}

gint glarea_key_press_event(GtkWidget *widget, GdkEventKey *event){
	switch (event->keyval) {
		case GDK_KEY_r:
			g_print("Button r pressed...redrawing\n");
			//gtk_widget_draw(glarea, (GdkRectangle *)NULL);
			break;
		case GDK_KEY_l:
			g_print("Button l pressed...redrawing\n");
			//gtk_widget_draw(glarea, (GdkRectangle *)NULL);
			break;
		case GDK_KEY_p:
			g_print("Button p pressed...redrawing\n");
			//gtk_widget_draw(glarea, (GdkRectangle *)NULL);
			break;
		case GDK_KEY_space:
		printf("WE HIT THE SPACE BAR\n");
			g_render_targets = !g_render_targets;
			break;
	}
	return (TRUE);
}

gboolean glarea_scroll_event (GtkWidget *widget, GdkEvent *event, gpointer data)
{
	//GdkDevice *source_device;
	GdkEventScroll *scroll_event;
	
	//source_device = gdk_event_get_source_device (event);
	scroll_event = (GdkEventScroll *) event;

	switch(scroll_event->direction){
		case GDK_SCROLL_UP: 
			break; 
		case GDK_SCROLL_DOWN: 
			break; 
    case GDK_SCROLL_LEFT:
    case GDK_SCROLL_RIGHT:
    case GDK_SCROLL_SMOOTH:
      break; 
	}
  return FALSE;
}

static GtkWidget *mk_spinner(const char *txt, GtkWidget *container,
                             float start, float min, float max, float step,
                             GtkCallback cb, gpointer data){
	GtkWidget *spinner, *label;
	GtkAdjustment *adj;
	GtkWidget *bx = gtk_box_new (GTK_ORIENTATION_HORIZONTAL, 1);

	label = gtk_label_new (txt);
	gtk_box_pack_start (GTK_BOX (bx), label, FALSE, TRUE, 2);
	gtk_widget_show(label);
	adj = (GtkAdjustment *)gtk_adjustment_new(
	          start, min, max, step, step, 0.0);
	float climb = 0.0;
	int digits = 0;
	if (step <= 0.001) {
		climb = 0.0001;
		digits = 4;
	} else if (step <= 0.01) {
		climb = 0.001;
		digits = 3;
	} else if (step <= 0.1) {
		climb = 0.01;
		digits = 2;
	} else if (step <= 0.99) {
		climb = 0.1;
		digits = 1;
	}
	spinner = gtk_spin_button_new (adj, climb, digits);
	gtk_spin_button_set_wrap (GTK_SPIN_BUTTON (spinner), FALSE);
	gtk_box_pack_start (GTK_BOX (bx), spinner, FALSE, TRUE, 2);
	g_signal_connect(spinner, "value-changed", G_CALLBACK(cb), data);
	gtk_widget_show(spinner);

	gtk_box_pack_start (GTK_BOX (container), bx, FALSE, FALSE, 2);
	return spinner;
}

GtkWidget *create_gl_window(){
	GtkWidget	*glwindow;

	/* ... Create the window which will hold glarea ... */
	/* ... Create a window in gtk - note the window is NOT visible yet ... */
	glwindow = gtk_window_new(GTK_WINDOW_TOPLEVEL);
	gtk_window_set_default_size (GTK_WINDOW (glwindow), 2500, 2100);

	/* ... Set window attributes ... */
	gtk_window_set_title(GTK_WINDOW(glwindow), "SHWFS");

	g_signal_connect_swapped (glwindow, "destroy",
							G_CALLBACK(destroyGUI), NULL);

	glarea = gtk_gl_area_new(); 
	
	gtk_gl_area_set_has_alpha((GtkGLArea*)glarea, TRUE); 
	gtk_gl_area_set_has_depth_buffer((GtkGLArea*)glarea, TRUE); 
	gtk_gl_area_set_auto_render((GtkGLArea*)glarea, TRUE); 
	g_signal_connect(glarea, "render", G_CALLBACK(render), NULL); 
	g_signal_connect(glarea, "realize", G_CALLBACK(on_realize), NULL); 
	if(!glarea) { 
		g_print("Can't create GtkGlArea widget\n");
		return FALSE;
	}
	
	GtkWidget *paned = gtk_paned_new (GTK_ORIENTATION_HORIZONTAL);
	gtk_container_add (GTK_CONTAINER (glwindow), paned);
	gtk_paned_add1(GTK_PANED (paned), glarea);

	/* ... Set up events and signals for OpenGL widget ... */
	gtk_widget_set_events(GTK_WIDGET(glarea),
								 gtk_widget_get_events (glarea) | 
								GDK_EXPOSURE_MASK|
								GDK_STRUCTURE_MASK | //resize? 
								GDK_BUTTON_PRESS_MASK|
								GDK_BUTTON_RELEASE_MASK|
								GDK_KEY_PRESS_MASK|
								GDK_KEY_RELEASE_MASK|
								GDK_POINTER_MOTION_MASK|
								GDK_POINTER_MOTION_HINT_MASK |
								GDK_SCROLL_MASK);
	g_signal_connect (glarea, "configure_event",
	                      G_CALLBACK(on_configure), NULL);
// 	g_signal_connect (glarea, "motion_notify_event",
// 								G_CALLBACK(glarea_motion_notify), NULL);
// 	g_signal_connect (glarea, "button_press_event",
// 								G_CALLBACK(glarea_button_press), NULL);
// 	g_signal_connect (glarea, "key_press_event",
// 								G_CALLBACK(glarea_key_press_event), NULL);
// 	g_signal_connect (glarea, "scroll-event",
// 								G_CALLBACK(glarea_scroll_event), NULL);

	//add in other widgets.
	GtkWidget *bx2 = gtk_box_new (GTK_ORIENTATION_VERTICAL, 1);
	
	g_centroidCalc_label.make("centroid calc time:", bx2); 
	g_framerate_label.make("framerate:", bx2); 
	g_dataSize_label.make("data size:", bx2); 
	
	GtkWidget *button = gtk_button_new_with_label("recalibrate");
	g_signal_connect(button,"clicked", 
							G_CALLBACK(recalibrate_clicked), NULL);
	gtk_box_pack_start (GTK_BOX (bx2), button, FALSE, FALSE, 2);
	
	button = gtk_button_new_with_label("reset recording");
	g_signal_connect(button,"clicked", 
							G_CALLBACK(reset_clicked), NULL);
	gtk_box_pack_start (GTK_BOX (bx2), button, FALSE, FALSE, 2);
	
	button = gtk_button_new_with_label("write data");
	g_signal_connect(button,"clicked", 
							G_CALLBACK(write_data_clicked), NULL);
	gtk_box_pack_start (GTK_BOX (bx2), button, FALSE, FALSE, 2);
	
	button = gtk_check_button_new_with_label("record data"); 
	g_signal_connect(button,"clicked", 
							G_CALLBACK(record_data_clicked), NULL);
	gtk_box_pack_start (GTK_BOX (bx2), button, FALSE, FALSE, 2);
	
	button = gtk_check_button_new_with_label("display centroids"); 
	gtk_toggle_button_set_active(
		GTK_TOGGLE_BUTTON(button), g_display_centroids); 
	g_signal_connect(button,"clicked", 
							G_CALLBACK(display_centroids_clicked), NULL);
	gtk_box_pack_start (GTK_BOX (bx2), button, FALSE, FALSE, 2);
	
	button = gtk_check_button_new_with_label("display camera image"); 
	gtk_toggle_button_set_active(
		GTK_TOGGLE_BUTTON(button), g_display_camimage); 
	g_signal_connect(button,"clicked", 
							G_CALLBACK(display_camimage_clicked), NULL);
	gtk_box_pack_start (GTK_BOX (bx2), button, FALSE, FALSE, 2);
	
	int i=0; 
	char lbl[128]; 
	snprintf(lbl, 128, "Exposure %d", i); 
	g_exposureSpin[i] = mk_spinner(lbl, bx2, 50, 50 , 1000, 10	, 
											updateExposureCB, GINT_TO_POINTER(0)); 
	
	
	g_logarea = gtk_text_view_new();
	GtkWidget* scrolledwindow = gtk_scrolled_window_new(NULL, NULL);
	gtk_container_add(GTK_CONTAINER(scrolledwindow), g_logarea);
	gtk_box_pack_start (GTK_BOX (bx2), scrolledwindow, TRUE, TRUE, 2);
	
	gtk_paned_add2 (GTK_PANED(paned), bx2);
	gint maxpain; 
	g_object_get(paned, "max-position", &maxpain, NULL); 
	// printf("max-position %d\n", maxpain); -- doesn't work, ?
	gtk_paned_set_position(GTK_PANED(paned), 2100); 

	gtk_widget_show_all(glwindow);
	/* ... Set focus to glarea widget and initialize OpenGL ... */
	gtk_widget_set_can_focus (glarea, TRUE); 
	gtk_widget_grab_focus(GTK_WIDGET(glarea));

	return(glwindow);
}

int main(int argc, char **argv){
	GtkWidget	*glwindow;
	
	mat_identity(g_mvp);  
	g_startTime = gettime(); 
	g_w = g_h = 2048; 
	for(int i=0; i<1; i++){
		g_s[i] = g_w * g_h; 
		g_data[i] = (unsigned char*)malloc(g_s[i]); 
		for(unsigned int j=0; j<g_s[i]/4; j++)
				g_data[i][j] = (unsigned char)(rand() & 0xff); 
		g_copy[i] = 1; 
	}

	
	/* ... Initialize gtk, handle command-line parameters ... */
	gtk_init(&argc, &argv);

	if(!(glwindow = create_gl_window())) {
		g_print("Can't create GtkGlArea widget/window!\n");
		return(0);
	}

	pthread_t thread1;
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	long i = 0; 
	i=0; pthread_create( &thread1, &attr, video_thread, (void*)i );
	
	// add a timeout to render frames, 60fps.  
	g_timeout_add (1000 / 25, redraw, 0);
	
	/* ... This is the event loop in gtk ... */
	/* ... Do not return until _gtk_main_quit_ is called ... */
	gtk_main();
	gl_fini(); 
	return(0);
}
