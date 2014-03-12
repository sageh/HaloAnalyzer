#include "common.h"

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <sstream>
#include <string>

#include <unistd.h>

// Helper function to check OpenGL attributes
void check_attribute(char *attname, int should_be, int attr) {
	if (attr != should_be) {
		fprintf(stderr, "Expected %s to be %d, found %d\n",
			attname, should_be, attr);
		exit(EXIT_FAILURE);
	}
}

// Get command line arguments pertaining to graphics.
int get_graphics_options(int argc, char **argv, int *oi, 
		int *w, int *h, int *bpp) {
	opterr = 0;
	// Use the given values as a starting point
	int wval = *w;
	int hval = *h;
	int dval = *bpp;
	int c;

	while ((c = getopt(argc, argv, "w:h:d:")) != -1) {
		switch (c) {
			case 'w':
				sscanf(optarg, "%d", &wval);
				break;
			case 'h':
				sscanf(optarg, "%d", &hval);
				break;
			case 'd':
				sscanf(optarg, "%d", &dval);
				break;
			case '?':
				if (isprint (optopt))
					fprintf(stderr, 
					"Unknown option '-%c'.\n", optopt);
				else
					fprintf(stderr,
					"Unknown option char '\\x%x'.\n",
					optopt);
				return 1;
			default:
				return 1;
		}
	}

	// Print what we got and set the values
	printf("getopt()'ed w = %d, h = %d, bpp = %d\n", wval, hval, dval);
	*w = wval;
	*h = hval;
	*bpp = dval;

	// Store the option index so the main program will know where
	// to continue.
	*oi = optind;
	return 0;
}

SDL_Surface* init_opengl_screen(int screen_width, int screen_height,
		int color_depth, int color_bits, int alpha_bits,
		int depth_bits, bool use_doublebuffering) {

	// The graphics context
	SDL_Surface *s;

	// Fetch video information (SDL_Init should already have been called)
	printf("Fetching SDL video information...\n");
	const SDL_VideoInfo *info = NULL;
	info = SDL_GetVideoInfo();

	if (!info) {
		fprintf(stderr, "Video information query failed: %s\n",
				SDL_GetError());
		exit(EXIT_FAILURE);
	}

	printf("SDL_GetVideoInfo suggests color depth of %d bpp\n",
		info->vfmt->BitsPerPixel);

	// Get all available video modes and print them
	printf("Fetching videomodes for SDL_OPENGL screen at %d bpp...\n",
		info->vfmt->BitsPerPixel);
	SDL_Rect **modes;
	modes = SDL_ListModes(info->vfmt, SDL_OPENGL);

	if (modes == NULL) {
		printf("No possible video modes!\n");
		exit(EXIT_FAILURE);
	}

	if (modes == (SDL_Rect**)-1) {
		printf("All resolutions available.\n");
	}
	else {
		printf("Found following modes:\n");
		for (int i=0; modes[i]; i++) 
			printf("%d x %d\n", modes[i]->w, modes[i]->h);
	}

	// Try to be very safe, by using the color depth provided by the
	// info struct if it differs from the one given as the argument
	if (color_depth != info->vfmt->BitsPerPixel) {
		printf( "Requested %d bpp, using SDL_GetInfo "
			"suggestion %d bpp\n", 
			color_depth, info->vfmt->BitsPerPixel);
		color_depth = info->vfmt->BitsPerPixel;
	}

	// Set up OpenGL attributes before creating the drawing context
	int doublebuffer = 0;
	if (use_doublebuffering)
		doublebuffer = 1;

	int result;
	printf( "Setting GL attributes:\n"
		"SDL_GL_RED_SIZE, %d\n"
		"SDL_GL_GREEN_SIZE, %d\n"
		"SDL_GL_BLUE_SIZE, %d\n"
		"SDL_GL_ALPHA_SIZE, %d\n"
		"SDL_GL_DEPTH_SIZE, %d\n"
		"SDL_GL_DOUBLEBUFFER, %d\n",
		color_bits, color_bits, color_bits, color_bits, depth_bits,
		doublebuffer);

	SDL_GL_SetAttribute(SDL_GL_RED_SIZE, color_bits);

	SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, color_bits);
	SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, color_bits);
	SDL_GL_SetAttribute(SDL_GL_ALPHA_SIZE, color_bits);
	SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, depth_bits);
	SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, doublebuffer);
	
	printf("Calling SDL_SetVideoMode(%d, %d, %d, SDL_OPENGL)\n",
		screen_width, screen_height, color_depth);
	s = SDL_SetVideoMode(screen_width, screen_height, color_depth, 
			SDL_OPENGL);

	// Check results
	SDL_GL_GetAttribute(SDL_GL_RED_SIZE, &result);
	check_attribute("SDL_GL_RED_SIZE", color_bits, result);

	SDL_GL_GetAttribute(SDL_GL_GREEN_SIZE, &result);
	check_attribute("SDL_GL_GREEN_SIZE", color_bits, result);

	SDL_GL_GetAttribute(SDL_GL_BLUE_SIZE, &result);
	check_attribute("SDL_GL_BLUE_SIZE", color_bits, result);
	
	SDL_GL_GetAttribute(SDL_GL_ALPHA_SIZE, &result);
	check_attribute("SDL_GL_ALPHA_SIZE", alpha_bits, result);

	SDL_GL_GetAttribute(SDL_GL_DEPTH_SIZE, &result);
	check_attribute("SDL_GL_DEPTH_SIZE", depth_bits, result);

	SDL_GL_GetAttribute(SDL_GL_DOUBLEBUFFER, &result);
	check_attribute("SDL_GL_DOUBLEBUFFER_SIZE", doublebuffer, result);

	// Check if we got a working context.
	if (s == NULL) {
		fprintf(stderr, "Unable to set video mode: %s\n",
				SDL_GetError());
		exit(EXIT_FAILURE);
	}
	printf("Set mode %d x %d at %d bpp\n", screen_width, screen_height,
			s->format->BitsPerPixel);

	// Initialize rest of OpenGL
	set_screen_size(screen_width, screen_height);

	// Set clearing color and clear screen
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT);

	// Set depth check
	glClearDepth(1.0);
	glDepthFunc(GL_LESS);
	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);

	// Be as accurate as possible
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);

	// Enable antialiasing
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);

	return s;
}

// Shape the GL window to size
void set_screen_size(int width, int height)
{
	glViewport(0,0,(GLsizei)(width),(GLsizei)(height));
}

void set_fov(double fov, int screen_width, int screen_height) {
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(fov,
		(GLfloat)(screen_width)/(GLfloat)(screen_height), 
		0.1f, 1000.0f);

	glMatrixMode(GL_MODELVIEW);
}

void set_text_mode(int screen_width, int screen_height) {
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(0, screen_width, 0, screen_height);
	glScalef(1, -1, 1);
	glTranslatef(0, -screen_height, 0);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
}

void unset_text_mode() {
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
}
	
void draw_box(double o[], double d[], GLfloat color[], 
		GLenum depthfunc) {

	// Store the old depth function, and activate the new one
	GLint old;
	glGetIntegerv(GL_DEPTH_FUNC, &old);
	glDepthFunc(depthfunc);

	// Set the color
	glColor4fv(color);

	// Draw the box as lines
	glBegin(GL_LINES);
	glVertex3d(0.0, 0.0, 0.0);
	glVertex3d(d[0], 0.0, 0.0);
	glVertex3d(0.0, 0.0, 0.0);
	glVertex3d(0.0, d[1], 0.0);
	glVertex3d(0.0, 0.0, 0.0);
	glVertex3d(0.0, 0.0, d[2]);
	glVertex3d(0.0, 0.0, d[2]);
	glVertex3d(d[0], 0.0,d[2]);
	glVertex3d(0.0, 0.0, d[2]);
	glVertex3d(0.0, d[1], d[2]);

	glVertex3d(d[0], d[1], d[2]);
	glVertex3d(0.0, d[1], d[2]); 
	glVertex3d(d[0], d[1], d[2]);
	glVertex3d(d[0], 0.0, d[2]); 
	glVertex3d(d[0], d[1], d[2]);
	glVertex3d(d[0], d[1], 0.0);
	glVertex3d(d[0], d[1], 0.0);
	glVertex3d(0.0, d[1], 0.0);
	glVertex3d(d[0], d[1], 0.0);
	glVertex3d(d[0], 0.0, 0.0);

	glVertex3d(d[0], 0.0, 0.0);
	glVertex3d(d[0], 0.0, d[2]);
	glVertex3d(0.0, d[1], 0.0);
	glVertex3d(0.0, d[1], d[2]);
	glEnd();

	// Restore old depth test function
	glDepthFunc(old);
}


int parse_halodump_file(vector<halo> *v, ifstream *f) {
	string row;
	istringstream ss;
	int current_mainhalo_id = -1;
	double tmp;
	halo main, sub;

	while (!f->eof()) {
		// Use istringstream to get the data
		getline(*f, row);
		ss.str(row);

		// First two fields are the main and subhalo indexes
		ss >> main.index >> sub.index;

		// npart, nvpart
		ss >> tmp >> tmp; 

		// The coordinates of main halo
		ss >> main.x >> main.y >> main.z;

		// 3 unnecessary fields
		for (int i=0; i < 3; i++)
			ss >> tmp;

		// Main halo m_vir, r_vir
		ss >> main.m_vir;
		ss >> main.r_vir;

		// Then 24 fields we do nothing with
		for (int i=0; i < 24; i++)
			ss >> tmp;

		// Subhalo x, y, z
		ss >> sub.x >> sub.y >> sub.z; 

		// 3 unnecessary fields
		for (int i=0; i < 3; i++)
			ss >> tmp;

		// Subhalo m_vir, r_vir
		ss >> sub.m_vir;
		ss >> sub.r_vir;

		// Then 30 more ignored fields
		for (int i=0; i < 30; i++)
			ss >> tmp;

		// Main halo formation time
		ss >> main.z_form;

		// Then 32 more ignored fields
		for (int i=0; i < 32; i++)
			ss >> tmp;

		// Sub halo formation time
		ss >> sub.z_form;

		// DEBUG: print the data
		/*
		printf("main: %d %f %f %f %f %f sub: %d %f %f %f %f %f\n",
		main.index, main.x, main.y, main.z, main.z_form, main.m_vir,
		sub.index, sub.x, sub.y, sub.z, sub.z_form, sub.m_vir);
		*/

		// Set the subhalo identifier boolean
		main.is_subhalo = false;
		sub.is_subhalo = true;

		// Only add the main halo if it is a different one
		if (main.index != current_mainhalo_id) {
			v->push_back(main);
			current_mainhalo_id = main.index;
		}

		// Always add the subhalo
		v->push_back(sub);

		// Double check for end of file
		if (f->eof())
			break;
	}

	// DEBUG:
	/*
	getline(*f, row);
	ss.str(row);
	int i=0;
	while (!ss.eof()) {
		ss >> tmp;
		printf("%d: %f\n", i, tmp);
	       i++;
	}	       
	*/

	return 1;
}

int parse_halos_file(vector<halo> *v, ifstream *f) {
	string row;
	istringstream ss;
	int halo_index = -1;
	double tmp;
	halo main;

	while (!f->eof()) {
		// Use istringstream to get the data
		getline(*f, row);
		ss.str(row);

		// Discard comments
		if (row.find("#", 0) != string::npos)
			continue;

		// Halo index is by row number
		halo_index++;
		main.index = halo_index;

		// npart, nvpart
		ss >> tmp >> tmp; 

		// The coordinates 
		ss >> main.x >> main.y >> main.z;

		// 3 unnecessary fields
		for (int i=0; i < 3; i++)
			ss >> tmp;

		// Main halo m_vir, r_vir
		ss >> main.m_vir;
		ss >> main.r_vir;


		// We have no information about subhalo status or
		// formation time, so we must set some default values.
		main.is_subhalo = false;
		main.z_form = 0.0;

		v->push_back(main);

		// Double check for end of file
		if (f->eof())
			break;
	}

	// DEBUG:
	/*
	getline(*f, row);
	ss.str(row);
	int i=0;
	while (!ss.eof()) {
		ss >> tmp;
		printf("%d: %f\n", i, tmp);
	       i++;
	}	       
	*/

	return 1;
}

int save_screenshot() {
	static int running_index = 0;

	// Construct file name
	static char filename[64];
	sprintf(filename, "screenshot_%d.bmp", running_index);

	// Save to file
	int result;

	SDL_Surface *image;
	SDL_Surface *temp;
	int idx;
	SDL_Surface * screen = SDL_GetVideoSurface();
	int w = screen->w;
	int h = screen->h;
	image = SDL_CreateRGBSurface(SDL_SWSURFACE, w, h, 24, 
			0x0000FF, 0x00FF00, 0xFF0000, 0x000000);
	temp = SDL_CreateRGBSurface(SDL_SWSURFACE, w, h, 24, 
			0x0000FF, 0x00FF00, 0xFF0000, 0x000000);
	glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, image->pixels);
	for (idx = 0; idx < h; idx++)
	{
		memcpy((unsigned char *)(temp->pixels) + 3 * w * idx,
			(unsigned char *)(image->pixels) + 3 * w * (h - idx),
				3*w);
	}
	memcpy(image->pixels,temp->pixels,w * h * 3);
	result = SDL_SaveBMP(image, filename);
	SDL_FreeSurface(image);
	SDL_FreeSurface(temp);

	running_index++;
	return (result == 0) ? 1 : 0;
}

