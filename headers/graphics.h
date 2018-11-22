#ifndef GRAPHICS_H
#define GRAPHICS_H

#include <GL/glew.h>
#include <GL/freeglut.h>
//#include <cuda_runtime.h>
//#include <cuda_gl_interop.h>
#include <vector>
#include <string>

#include "../utils/simple_timer.h"
#include "../utils/simple_initializer.h"
#include "../utils/simple_palettes.h"

#include "../glm/glm.hpp"
#include "../glm/gtc/matrix_transform.hpp"
#include "../glm/gtc/type_ptr.hpp"

using namespace std;


class Palette{
	public:
	int n;
	vector <glm::vec4> colors;

	Palette(int n);
	void create_rainbow(float start=0, float end=0.75);
	void create_random(float start=0, float end=1);
	void create_grayscale(float start=0, float end=1);
	void create_ramp(glm::vec4 start, glm::vec4 end);	
	void print();
	vector <float> map_values(float * v, int nval, int stride = 1, float vmin = 1e20, float vmax = 1e20);
};

/* =======================================================================
	Shape class
======================================================================= */ 
class Shape{
	public:
	string objName;
	string type;
//	bool doubleBuffered;
	int nVertices;	//!< Number of vertices that comprise the shape
	int dim; 		//!< Number of components per vertex (2 for 2D, 3 for 3D)
	
	GLuint vbo_id;
	GLuint colorBuffer_id;

	GLuint vertexShader_id;
	GLuint fragmentShader_id;
	GLuint program_id;
	
	string vertexShaderFile;
	string fragmentShaderFile;
	
	glm::mat4 model;
	
	float pointSize;
	
	bool b_render;
	
	public:
//	Shape(){};
	Shape(int nVert, int components_per_vertex, string _type, bool ren = true);
	~Shape();

	void setVertices(void* data);
	void createShaders();
	void createColorBuffer();
	
	void autoExtent(float* data);
	void setExtent(vector <float>& ex);
	
	void deleteVBO();
	void deleteShaders();
	void deleteColorBuffer();
	
	void setColors(float* colData);
	
	void setRenderVariable(string s, float  f);
	void setRenderVariable(string s, glm::vec2 f);
	void setRenderVariable(string s, glm::vec3 f);
	void setRenderVariable(string s, glm::vec4 f);
	void setShaderVariable(string s, glm::mat4 f);

	void render();
	
	void useProgram();
};


class Shape2D : public Shape{
	public:
	Shape2D(int nVert, string _type);
	void setExtent(float xmin, float xmax, float ymin, float ymax);
};

//class ColorMap : public Shape{
//	public:
//	int nlevs;
//	int nx;
//	float xmin, xmax;

//	vector <Colour_rgb> palette;

//	ColorMap();
//	ColorMap(string obj_name, int nlevs, int nx, float _xmin, float _xmax);
//	
//	void createGridCentres(float2* cmap_tmp);
//	void updateColors(float* colData, int nCol, float cmin = -1e20, float cmax = -1e20);
//	void render();

//	void destroy();
//};



/* =======================================================================
	NOTE ON USING A RENDERER

	Renderer should be used as follows:

	Renderer R;
	R.init();	// renderer init MUST happen before GL init.
	initGL();	// requires valid initialized renderer to be set for use by openGL 
	
	R.connect(particle_system_1);
		... render
	R.disconect();
	
	R.connect(particle_system_2);
		... render
	R.disconnect();
	
	...

=======================================================================  */ 


/* =======================================================================
	NOTE ON PARTICLE COLURING METHOD
	
	class Particle has various particle attributes starting wA.
	All attributes are of type float or int (4 bytes).
	particleColorAttribute is a void pointer to the desired attribute 
		in the first particle. i.e. &pvec[0].<attr>
	since all particles are 4 bytes, this pointer can be made to point 
		at the nth attribute as (float*)&pvec[0].<attr>+n. This is how 
		the current coloring attributes are set using the number keys.
		the value 'n' is stored in the iColorAttribute
	once the pointer is set at the first value, it can be advanced by 
		sizeof(Particle) bytes to get the same attribute value in the 
		next particle. This is used to set the color buffer.
	Of course, we need to know the type of the attribute to do this, 
		so we store the types (int/float) in an array.
	To get the final color, the value is discretized using the min/max
		bounds and the palette is looked up to get RGB.
=======================================================================  */ 


enum UpdateMode {Step, Time};

class Renderer{
	private:

	// string to store command to be executed from GUI console
	string command;

	public:	
	
	glm::mat4 view, projection;
	float camera_tx, camera_ty;
	float camera_s, camera_rx, camera_ry;
	
	int up_axis;
	
	// update intervals/steps
	int nSkip;				// number of steps to skip before re-rendering
	int displayInterval;	// interval in ms between 2 display calls
	int updateMode; 		// update mode: update after fixed time or fixed steps
	int quality;			// quality of graphics
	bool b_paused;
	//int b_anim_on;

	// layers - by default only layer 0 is visible
	vector <bool> layerVis;
	
//	// colour palettes
//	vector <Colour_rgb> palette;
//	vector <Colour_rgb> palette_rand;
//	vector <Colour_rgb> palette_grad;

	unsigned int window_width;
	unsigned int window_height;

	// coordinate system
	float xmin, xmax, ymin, ymax, zmin, zmax;

	// render flags
	bool b_renderConsole;
	bool b_renderText;
	bool b_renderLabels;
	bool b_renderRs;
	bool b_renderGrid;
	bool b_renderAxes;
	bool b_renderColorMap;
		
	// frame counter
	SimpleCounter frameCounter;
	
	// GPU stuff
	GLuint vao_id;

	int swap;	// index of the most recently updated buffer	

	vector <Shape*> shapes_vec;	// shapes to render


	public:
	// init
	void init();

	// fancy stuff
	int getDisplayInterval();

	// add shapes to render list
	int addShape(Shape* shp);

	int  renderConsole();
	int  renderAxes(float lim, float trans);
	void renderGrid();
	
	void receiveConsoleChar(char key);
	int executeCommand();
	string makeTitle();
	
	void togglePause();
	void toggleConsole();
	void toggleText();
	void toggleGrid();
	void toggleAxes();
};

// pointers to the particle system to display and renderer to render display
extern Renderer * glRenderer;
extern int generic_count;
//extern ParticleSystem * glPsys;

// util functions
void loadShader(string filename, GLuint &shader_id, GLenum shader_type);
vector <float> calcExtent(float* data, int nVertices, int dim);

// openGL callbacks
bool init_hyperGL(int *argc, char **argv);
void timerEvent(int value);
void reshape(int w, int h);
void keyPress(unsigned char key, int x, int y);
void specialKeyPress(int key, int x, int y);
void mouseMove(int x, int y);
void mousePress(int button, int state, int x, int y);
void display();
void cleanup();


#endif


