/* Common functionality shared by all visualization programs.
 *
 */

#ifndef COMMON_H
#define COMMON_H

#include <GL/gl.h>

#include <iostream>
#include <fstream>
#include <vector>

#include <SDL.h>

using namespace std;


//
// Datatype definitions
///////////////////////

// A halo structure
typedef struct _halo {
	// Halo index
	int index;

	// Coordinates
	double x, y, z;

	// Formation time as a redshift
	double z_form;

	// Virial mass and radius
	double m_vir;
	double r_vir;

	// Whether this halo is a subhalo
	bool is_subhalo;
} halo;

/* Get the command line arguments that have to do with setting the graphics
 * mode, and return the index from where to continue. Returns 0 if all went
 * well, nonzero otherwise. */
int get_graphics_options(int argc, char **argv, int *index,
		int *w, int *h, int *bpp);

/* Initialize a SDL graphics context for OpenGL use. */
SDL_Surface* init_opengl_screen(int screen_width, int screen_height, 
		int color_depth, int color_bits, int alpha_bits, 
		int depth_bits, bool double_buffer);

/* (Re)set OpenGL screen size. */
void set_screen_size(int width, int height);

/* Set Field of View. */
void set_fov(double fov, int screen_width, int screen_height);

/* Functions for setting up a mode suitable for drawing text (but not much else)
 * and reversing it. */
void set_text_mode(int screen_width, int screen_height); 
void unset_text_mode();

/* Draw a box using the given origin, dimensions, color (Glfloat[4]) and
 * GL depth test function.
 * Implicitly assumes that origin and dimension contain 3 values and color 4.*/
void draw_box(double origin[], double dimensions[], GLfloat color[], 
		GLenum depth_function = GL_ALWAYS);

/* Parse a halodump_halos_with_subhalos.dat file that must have additional
 * formation time redshift and data columns added. Stores the data into
 * the given vector of halo structs. Returns true on success, false otherwise.*/
int parse_halodump_file(vector<halo> *v, ifstream *f);

/* Parse a _halos file into the given vector of halo structs. */
int parse_halos_file(vector<halo> *v, ifstream *f);

/* Grabs a screenshot and saves it to a file. Returns true on success, false
 * otherwise.
 * Uses sequential numbering for screenshots taken in one session. Overwrites
 * existing files. */
int save_screenshot();


#endif /* common.h */
