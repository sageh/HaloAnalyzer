// No point in not using opengl
#define USE_OPENGL	1

#ifdef USE_OPENGL
// OpenGL
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// SDL
#include <SDL.h>

#include "common.h"

using namespace std;

// Datatype defs
////////////////

typedef struct _camera {
	double x, y, z;
	double phi, theta;
	double r;
	double up_x, up_y, up_z;
} camera;

// Ugly globals
//////////////
camera cam;

// Halo data for each frame
vector<vector<halo>*> halo_frames;
// Redshift values for each frame
vector<string> redshifts;

double box_size;
double field_of_view = 40.0;
int mouse_mode = 0;
unsigned int current_frame = 0;

// Is the animation cycling
bool animate = false;

// Mouse motion y-axis invert
bool invert_mouse = false;

// Whether to print parameter values etc. 
bool display_info = true; 

// Constant defs
////////////////

// Screen and video mode specs
int screen_width = 1024;
int screen_height = 768;
int color_bits = 8;
int colordepth_bits = 32;
int depth_bits = 16;
int doublebuffer = 1;
const int fps_cap = 20;

// Input handling
const double mouse_constant = 100;
const double keymove_constant = 0.1;
const int MOUSE_ROTATE 	= 1;
const int MOUSE_ZOOM 	= 2;

// Graphic look
double z_bg_cutoff = 1.0;
double z_gr_cutoff = 2.0;
double z_r_factor = 0.01;
double z_g_factor = 0.01;
double z_b_factor = 0.2;
double halo_alpha = 0.7;
double mass_zero_point = 8.0;

// Animation
const int FRAME_TICKS = 1000;

// Misc.
const unsigned int max_strlen = 256;

// Function defs
////////////////
void usage();
void handle_key_event(SDL_Event *e);
void handle_mouse_click(SDL_Event *e);
void handle_mouse_motion(SDL_Event *e);
void calculate_camera(camera *c);
void set_mouse_mode(int mode, bool on_or_off);
void print_info();
void draw_halos(vector<halo> *halos);

// For sorting halos into descending order by distance from the camera
void sort_halos_by_depth(vector<halo> *h);
bool halo_depth_cmp(halo a, halo b) {
	return (a.x-cam.x)*(a.x-cam.x) 
		+ (a.y-cam.y)*(a.y-cam.y) 
		+ (a.z-cam.z)*(a.z-cam.z) >
		(b.x-cam.x)*(b.x-cam.x) 
		+ (b.y-cam.y)*(b.y-cam.y) 
		+ (b.z-cam.z)*(b.z-cam.z);
}


// For GLUT bitmap strings
void bitmap_string(float x, float y, void *font, char *string) {  
	char *c;
	glRasterPos2f(x, y);
	for (c=string; *c != '\0'; c++) {
		glutBitmapCharacter(font, *c);
	}
}

// Calculate a color based on the halo mass
void calculate_color(GLdouble *color, double m_vir) {
	// Use the 10-base logarithm as a basis, with variable zero point
	double log_m = log10(m_vir);
	log_m -= mass_zero_point;

	/* Whether to cycle from one color to another, or to adjust all at
	 * once. */
#ifdef ROTATING_COLORS
	// If result is negative, then set to blue
	if (log_m < 0.0) {
		color[0] = 0.0;
		color[1] = 0.0;
		color[2] = 1.0;
		return;
	}
	if (log_m < z_bg_cutoff) {
		color[0] = 0.0;
		color[1] = log_m/z_bg_cutoff;
		color[2] = 1.0;
	}
	else if (log_m < z_gr_cutoff) {
		color[0] = (log_m-z_bg_cutoff)/(z_gr_cutoff-z_bg_cutoff);
		color[1] = 1.0;
		color[2] = 1.0 - (log_m-z_bg_cutoff)/(z_gr_cutoff-z_bg_cutoff);
	}
	else {
		color[0] = 1.0; 
		color[1] = exp(-(log_m-z_gr_cutoff)/z_gr_cutoff);
		color[2] = 0.0;
	}
#else
	// If result is negative, then set to 0
	if (log_m < 0.0) {
		color[0] = color[1] = color[2] = 1.0;
		return;
	}

	color[0] = exp(-z_r_factor*log_m);
	color[1] = exp(-z_g_factor*log_m);
	color[2] = exp(-z_b_factor*log_m);
#endif
}

// Calculate size (in this case, radius) from halo mass
double calculate_size(halo *h) {
	return h->r_vir/1000.0;
}

void render(SDL_Surface *screen, vector<halo> *halos, camera *cam) {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// Calculate new camera parameters
	calculate_camera(cam);

	// Look at the center of the box (assume 10 MPc for now)
	gluLookAt(cam->x, cam->y, cam->z, 
			box_size/2, box_size/2, box_size/2,
			cam->up_x, cam->up_y, cam->up_z);

	// Draw the box
	double origin[] = {.0, .0, .0};	
	double dim[] = {box_size, box_size, box_size};
	static GLfloat boxcol[] = {1.0, 1.0, 1.0, 1.0};
	draw_box(origin, dim, boxcol);

	// Before drawing the halos, we should sort them.
	sort_halos_by_depth(halos);

	// Draw all halos
	draw_halos(halos);

	glFlush();

	// Set a mode that gets us the familiar screen coordinates
	set_text_mode(screen_width, screen_height);

	// Draw mass scale, first the colors, then texts.
	// For the scale, we need to have blending on, so that the colors
	// correspond to what is seen on the screen.
	GLfloat x, y;
	GLfloat w = 20;
	GLfloat h = 10;
	GLdouble color[3];
	x = 10;
	glEnable(GL_BLEND);
	glBegin(GL_QUADS);
	for (int i=0; i < 30; i++) {
		y = 200 + i*h;
		calculate_color(color, pow(10, i/5.0 + mass_zero_point));
		glColor4d(color[0], color[1], color[2], 1.0);
		glVertex3f(x, 	y, 	0.0f);
		glVertex3f(x+w, y, 	0.0f);
		glVertex3f(x+w,	y+h, 	0.0f);
		glVertex3f(x, 	y+h, 	0.0f);
	}
	glEnd();
	glDisable(GL_BLEND);

	// Set bright white for text
	glColor4d(1.0, 1.0, 1.0, 1.0);
	// String buffer
	char str[max_strlen];
	for (int i=0; i < 30; i++) {
		y = 200 + i*h;
		snprintf(str, max_strlen, "m_vir = %2.2e", 
				pow(10, i/5.0 + mass_zero_point));
		bitmap_string(x+w+2, y+h, GLUT_BITMAP_HELVETICA_10, str);
	}

	// Print the parameter values and other misc info.
	print_info();

	unset_text_mode();

	glFlush();
	SDL_GL_SwapBuffers();
}

int main(int argc, char **argv) {
	// Our drawing context
	SDL_Surface *screen;

	// Set camera to some sensible defaults
	cam.x = cam.y = cam.z = 0.0;
	cam.phi = 0.0;
	cam.theta = M_PI/2;
	cam.r = 10;
	cam.up_x = cam.up_y = 0.0;
	cam.up_z = 1.0;
	calculate_camera(&cam);

	// Check arguments and read the data

	// First get the graphics arguments.
	int oi = 0; 
	if (get_graphics_options(argc, argv, &oi, 
		&screen_width, &screen_height, &colordepth_bits) != 0) {
		// Error, bail out.
		printf("getopt() failed.\n");
		usage();
		exit(EXIT_FAILURE);
	}
	printf("Resuming parsing arguments from index %d\n", oi);

	// In addition to graph arguments, we also need some others. 
	if (argc - oi < 2) {
		printf("Not enough arguments.\n");
		usage();
		exit(EXIT_FAILURE);
	}

	box_size = strtod(argv[oi], NULL);
	if (box_size < 0.1) {
		usage();
		return 0;
	}

	printf("Box size is %lf megaparsecs\n", box_size);
	cam.r = box_size;

	// Read one halo data file for each frame
	for (int i=oi+1; i < argc; i++) {
		string fname(argv[i]);
		vector<halo> *input = new vector<halo>();

		// Read halo data
		printf("Opening input file %s\n", fname.c_str());
		ifstream fin(fname.c_str(), ios::in);
		if (!fin) {
			fprintf(stderr, "Error opening file: %s\n", 
					fname.c_str());
			delete input;
			return 1;
		}


		printf("Reading in halo data\n");
		if (!parse_halos_file(input, &fin)) {
			fprintf(stderr, "Error reading data from file: %s\n",
					fname.c_str());
			delete input;
			fin.close();
			return 1;
		}
		fin.close();
		printf("Read in data for %d halos\n", input->size());

		// The read was a success so strip the redshift value and
		// store it. Assume that the redshift value is exactly 5
		// characters long.
		unsigned int z_start = fname.find("SIMU.z", 0);
		if (z_start == string::npos) {
			fprintf(stderr, "Couldn't find redshift value from "
					"file name: %s\n", fname.c_str()); 
			delete input;
			return 1;
		}
		string redshift = fname.substr(z_start+6, 5);

		// Store the halo data and the redshift value
		halo_frames.push_back(input);
		redshifts.push_back(redshift);
	}

	// Initialize SDL library
	/////////////////////////
	printf("SDL_Init\n");
	if (SDL_Init(SDL_INIT_EVERYTHING) < 0) {
		fprintf(stderr, "Unable to init SDL: %s\n", SDL_GetError());
		exit(1);
	}
	atexit(SDL_Quit);

	// Set keyboard repeat
	SDL_EnableKeyRepeat(SDL_DEFAULT_REPEAT_DELAY, 
			SDL_DEFAULT_REPEAT_INTERVAL);

	// Initialize our drawing context
	printf("Initializing screen\n");
	screen = init_opengl_screen(screen_width, screen_height,
			colordepth_bits, color_bits, color_bits, 
			depth_bits, true);
	set_fov(field_of_view, screen_width, screen_height);

	// Loop until user quits 
	bool quitflag = false;
	int start_ticks;
	long int last_frame_ticks = 0;
	while (!quitflag)
	{
		// Start measuring the time it takes to draw this frame 
		start_ticks = SDL_GetTicks();

		// If we are animating, and it is time for a new frame, then
		// switch the data.
		if (animate && start_ticks > last_frame_ticks) {
			last_frame_ticks = start_ticks + FRAME_TICKS;
			if (++current_frame >= halo_frames.size())
				current_frame = 0;

			/*
			printf("frame %d of %d\n", current_frame+1, 
					halo_frames.size());
					*/
		}

		// Poll for events, and handle the ones we care about.
		SDL_Event event;
		while (SDL_PollEvent(&event)) 
		{
			switch (event.type) 
			{
				case SDL_KEYDOWN:
					handle_key_event(&event);
					break;
				case SDL_MOUSEBUTTONDOWN:
				case SDL_MOUSEBUTTONUP:
					handle_mouse_click(&event);
					break;
				case SDL_MOUSEMOTION:
					handle_mouse_motion(&event);
					break;
				case SDL_QUIT:
					quitflag = true;
					break;
			}
		}



		// Render screen
		//render(screen, &halos, &cam);
		render(screen, halo_frames[current_frame], &cam);

		// If there is extra time, just wait
		while ((SDL_GetTicks() - start_ticks) < 1000 / fps_cap);
	}

	// Possible cleanup here

	return 0;
}

void usage() {
	printf( "Arguments: GRAPHICS_OPTS box size (in MPc), "
		"halos files (descending z)\n"
		"GRAPHICS_OPTS:\n"
		"-w screen width\n"
		"-h screen height\n"
		"-d color depth in bits\n"
		);
}

void keep_positive(double *d) {
	if (*d <= 0.0)
		*d = FLT_EPSILON;
}

void keep_non_negative(double *d) {
	if (*d < 0.0)
		*d = 0;
}

void keep_within(double *d, double low, double high) {
	if (*d < low)
		*d = low;
	else if (*d > high)
		*d = high;
}

void handle_key_event(SDL_Event *event) {
	// If escape is pressed, return (and
	// thus, quit)
	static SDL_Event quit;
	static const double color_delta = 0.005;
	static const double alpha_delta = 0.05;
	static const double fov_delta 	= 1.00;
	static const double cutoff_delta = 0.25;
	static const double mass_delta = 0.25;

	// Get the modifider key status
	SDLMod mods = event->key.keysym.mod;
	switch (event->key.keysym.sym) {
		case SDLK_ESCAPE:
			// Queue quit
			quit.type = SDL_QUIT;
			SDL_PushEvent(&quit);
			break;
		case SDLK_SPACE:
			animate = !animate;
			break;
		case SDLK_LEFT:
			if (current_frame-- == 0)
				current_frame = halo_frames.size()-1;
			break;
		case SDLK_RIGHT:
			if (++current_frame >= halo_frames.size())
				current_frame = 0;
			break;
		case SDLK_UP:
			cam.r -= keymove_constant;
			break;
		case SDLK_DOWN:
			cam.r += keymove_constant;
			break;
		case SDLK_r:
			if (mods & KMOD_SHIFT)
				z_r_factor -= color_delta;
			else
				z_r_factor += color_delta;
			keep_non_negative(&z_r_factor);
			break;
		case SDLK_g:
			if (mods & KMOD_SHIFT)
				z_g_factor -= color_delta;
			else
				z_g_factor += color_delta;
			keep_non_negative(&z_g_factor);
			break;
		case SDLK_b:
			if (mods & KMOD_SHIFT)
				z_b_factor -= color_delta;
			else
				z_b_factor += color_delta;
			keep_non_negative(&z_b_factor);
			break;
		case SDLK_m:
			if (mods & KMOD_SHIFT)
				mass_zero_point -= mass_delta;
			else
				mass_zero_point += mass_delta;
			keep_non_negative(&z_b_factor);
			break;
		case SDLK_a:
			if (mods & KMOD_SHIFT)
				halo_alpha -= alpha_delta;
			else
				halo_alpha += alpha_delta;
			keep_within(&halo_alpha, 0.0, 1.0);
			break;
		case SDLK_f:
			if (mods & KMOD_SHIFT)
				field_of_view -= fov_delta;
			else
				field_of_view += fov_delta;
			keep_within(&field_of_view, 1.0, 120.0);
			set_fov(field_of_view, screen_width, screen_height);
			break;
		case SDLK_1:
			if (mods & KMOD_SHIFT)
				z_bg_cutoff -= cutoff_delta;
			else
				z_bg_cutoff += cutoff_delta;
			keep_within(&z_bg_cutoff, 0.0, z_gr_cutoff);
			break;
		case SDLK_2:
			if (mods & KMOD_SHIFT)
				z_gr_cutoff -= cutoff_delta;
			else
				z_gr_cutoff += cutoff_delta;
			keep_within(&z_bg_cutoff, z_bg_cutoff, 100.0);
			break;
		case SDLK_i:
			invert_mouse = !invert_mouse;
			break;
		case SDLK_d:
			display_info = !display_info;
			break;
		case SDLK_F12:
			save_screenshot();
			break;
		default:
			break;
	}
}

void set_mouse_mode(int mode, bool turn_on) {
	if (turn_on) {
		// Grab input and hide the mouse cursor
		SDL_ShowCursor(0);
		SDL_WM_GrabInput(SDL_GRAB_ON);
		mouse_mode |= mode;
	}
	else {
		// Unhide cursor, and yield input
		// Grab input and hide the mouse cursor
		SDL_ShowCursor(1);
		SDL_WM_GrabInput(SDL_GRAB_OFF);
		mouse_mode &= ~mode;
	}
}

void handle_mouse_click(SDL_Event *event) {
	SDL_MouseButtonEvent m = event->button;
	switch (m.button) {
		case SDL_BUTTON_LEFT:
			if (m.type == SDL_MOUSEBUTTONDOWN)
				set_mouse_mode(MOUSE_ROTATE, true);
			else
				set_mouse_mode(MOUSE_ROTATE, false);
			break;
		case SDL_BUTTON_RIGHT:
			if (m.type == SDL_MOUSEBUTTONDOWN)
				set_mouse_mode(MOUSE_ZOOM, true);
			else
				set_mouse_mode(MOUSE_ZOOM, false);
			break;
	}
}

void handle_mouse_motion(SDL_Event *event) {
	// If a button is pressed, then zoom, otherwise rotate
	/*
	SDL_PumpEvents();
	if (SDL_GetMouseState(NULL, NULL) & SDL_BUTTON(1)) {
		cam.r += event->motion.yrel/mouse_constant;
		keep_positive(&cam.r);
	}
	else {
		// Rotate camera
		cam.phi += event->motion.xrel/mouse_constant;
		cam.theta += event->motion.yrel/mouse_constant;
	}
	*/
	if (mouse_mode & MOUSE_ROTATE) {
		// Rotate camera
		cam.phi += event->motion.xrel/mouse_constant;
		if (invert_mouse)
			cam.theta -= event->motion.yrel/mouse_constant;
		else
			cam.theta += event->motion.yrel/mouse_constant;
	}

	if (mouse_mode & MOUSE_ZOOM) {
		cam.r += event->motion.yrel*box_size/(10.0*mouse_constant);
		keep_positive(&cam.r);
	}
}

void calculate_camera(camera *cam) {
	// Calculate new position
	cam->x = box_size/2 + cam->r*sin(cam->theta)*cos(cam->phi);
	cam->y = box_size/2 + cam->r*sin(cam->theta)*sin(cam->phi);
	cam->z = box_size/2 + cam->r*cos(cam->theta);

	// Calculate where the up direction is
	cam->up_x = -cos(cam->theta)*cos(cam->phi);
	cam->up_y = -cos(cam->theta)*sin(cam->phi);
	cam->up_z = sin(cam->theta);
}

void sort_halos_by_depth(vector<halo> *halos) {
	sort(halos->begin(), halos->end(), halo_depth_cmp);
}

void draw_halos(vector<halo> *halos) {
	// Draw every halo
	// If we don't have translucency, then don't blend and enable depth
	// testing
	if (halo_alpha < 1.0) {
		glEnable(GL_BLEND);
		glDepthFunc(GL_ALWAYS);
	}
	else {
		glDisable(GL_BLEND);
		glDepthFunc(GL_LESS);
	}

	// Use source alpha for blending, because destination alpha blending
	// doesn't work with at least M$:s software opengl renderer
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	//glBlendFunc(GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA);
	GLUquadricObj *quadric = gluNewQuadric();
	GLdouble color[3];
	GLdouble size;
	for (unsigned int i=0; i < halos->size(); i++) {
		halo h;
		h = (*halos)[i];

		// Calculate point color from the halo formation age
		calculate_color(color, h.m_vir);

		// Calculate size from the mass
		size = calculate_size(&h);;

		// Draw a gluSphere
		//glColor3dv(color);
		glColor4d(color[0], color[1], color[2], halo_alpha);
		glPushMatrix();
		glTranslated(h.x, h.y, h.z);
		gluSphere(quadric, size, 8, 6);
		glPopMatrix();

	}
	glDisable(GL_BLEND);
}

void print_info() {
	// Don't print, if we want an uncluttered screen
	if (!display_info)
		return;

	// Use a static buffer
	static char str[max_strlen];

	// Print the current redshift
	snprintf(str, max_strlen, "Redshift: %s", 
			redshifts[current_frame].c_str());
	bitmap_string(screen_width/2-100, 11, GLUT_BITMAP_HELVETICA_10, str);

	// Print r, g, b and alpha color parameter values
	snprintf(str, max_strlen, "Red: %f", z_r_factor);
	bitmap_string(0, 11, GLUT_BITMAP_HELVETICA_10, str);
	snprintf(str, max_strlen, "Green: %f", z_g_factor);
	bitmap_string(0, 22, GLUT_BITMAP_HELVETICA_10, str);
	snprintf(str, max_strlen, "Blue: %f", z_b_factor);
	bitmap_string(0, 33, GLUT_BITMAP_HELVETICA_10, str);
	snprintf(str, max_strlen, "Alpha: %f", halo_alpha);
	bitmap_string(0, 44, GLUT_BITMAP_HELVETICA_10, str);

	// Print blue-green, green-red cutoff values
	snprintf(str, max_strlen, "BG cutoff: %f", z_bg_cutoff);
	bitmap_string(0, 55, GLUT_BITMAP_HELVETICA_10, str);
	snprintf(str, max_strlen, "GR cutoff: %f", z_gr_cutoff);
	bitmap_string(0, 66, GLUT_BITMAP_HELVETICA_10, str);

	// Print FOV etc.
	snprintf(str, max_strlen, "FOV: %f", field_of_view);
	bitmap_string(0, 77, GLUT_BITMAP_HELVETICA_10, str);


	// XXX: Print some debug info
	snprintf(str, max_strlen, "mouse_mode: %d", mouse_mode);
	bitmap_string(0, 88, GLUT_BITMAP_HELVETICA_10, str);
	snprintf(str, max_strlen, "invert_mouse: %d", invert_mouse);
	bitmap_string(0, 99, GLUT_BITMAP_HELVETICA_10, str);
	snprintf(str, max_strlen, "animate: %s", animate ? "true" : "false");
	bitmap_string(0, 110, GLUT_BITMAP_HELVETICA_10, str);
}


