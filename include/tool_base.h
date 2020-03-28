#ifndef QUICKGL_TOOL_BASE_H_
#define QUICKGL_TOOL_BASE_H_

#include <list>
#include <iostream>

#include <glm/glm.hpp>
//#include "../include/shape.h"
#include  "../include/glinit.h"

//class Tool{
//	public:
//	Shape * affectedShape; 		// affected shape
//	Shape * displayShape; 		// tool's visualization shape

//	public:
//	Tool(Shape * _affectedShape);
//	virtual ~Tool(){};
//	virtual void update();

//	virtual int  cursorLocation(int x, int y){}; // tells if cursor is outside (0), inside(1), bottom-edge (21), left-edge (22), top edge(23), right edge (24)
//	virtual void changeCursor(int x, int y){};

//	virtual bool blockNewSelection(int x, int y){};
//	
//	virtual void initialize(int x, int y){};
//	virtual void transform(int x, int y){};
//	virtual void finalize(int x, int y){};
//	virtual void mouseMove(int x, int y);
//	virtual void mouseHover(int x, int y);
//	virtual void mouseClick(int x, int y){};
//};



// QuickGL follows the principle of creating all objects (including cameras) aligned to axes and then using the matrices to get then in position

class Tool{
	protected:
	bool lClick, mClick, rClick;	// whether any mouse button was clicked
	bool wheelUp, wheelDown;		// whether wheel was moved (user is responsible to reset this to false)
	int x0, y0;						// mouse coordinates when any button was clicked

	public:
	static std::list <Tool*> activeTools;
	
	public:
	Tool();	
	virtual ~Tool();

	void activate();
	void captureClick(int button, int state, int mods);	

	virtual void onClick(int button, int state, int mods);
	virtual void onMouseMove(int x, int y);
	virtual void onScroll(double dx, double dy);
};


void onClick1(GLFWwindow* window, int button, int state, int mods);

void onMouseMove1(GLFWwindow* window, double x, double y);

void onScroll1(GLFWwindow* window, double dx, double dy);






#endif
