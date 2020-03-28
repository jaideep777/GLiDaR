#include  "../include/tool_base.h"

#include <algorithm>

//#include "../include/camera.h"
using namespace std;


list <Tool*> Tool::activeTools;


Tool::Tool(){
	mClick = lClick = rClick = 0;
	activeTools.push_front(this);
	cout << "Constructing tool #" << activeTools.size() << ", " << activeTools.size() << " in list" << endl;

	// this is a bit suboptimal, but including here to prevent tool dependency in glinit
    glfwSetMouseButtonCallback(window, onClick1);
	glfwSetCursorPosCallback(window, onMouseMove1);
	glfwSetScrollCallback(window, onScroll1);

}


Tool::~Tool(){
	try{
		list<Tool*>::iterator it = find(activeTools.begin(), activeTools.end(), this);
		if (it == activeTools.end()) throw string("Tool destructor could not find tool in active list!");
		activeTools.erase(it);
		cout << "Destroying tool, " << activeTools.size() << " left in list." << endl;
	}
	catch(string s){ 
		cout << "FATAL ERROR: " << s << endl; 
	}
}


void Tool::activate(){
	// bring to front of queue
	try{
		list<Tool*>::iterator it = find(activeTools.begin(), activeTools.end(), this);
		if (it == activeTools.end()) throw string("Tool destructor could not find tool in active list!");
		activeTools.splice(activeTools.begin(), activeTools, it);	// erase *it from current location (it) and put it at begin()
	}
	catch(string s){ 
		cout << "FATAL ERROR: " << s << endl; 
	}
}


void Tool::captureClick(int button, int state, int mods){
	double x = 0, y=0;
	glfwGetCursorPos(window, &x, &y);

	if (button == GLFW_MOUSE_BUTTON_LEFT){
		if (state == GLFW_PRESS){
			cout << "L Click captured by Tool at: " << x << " " << y << endl;
			lClick = true;
			x0 = x;
			y0 = y;
		}
		else{
			lClick = false;
		} 
	}

	if (button == GLFW_MOUSE_BUTTON_MIDDLE){
		if (state == GLFW_PRESS){
			cout << "M Click captured by Tool at: " << x << " " << y << endl;
			mClick = true;
			x0 = x;
			y0 = y;
		}
		else{
			mClick = false;
		} 
	}

	if (button == GLFW_MOUSE_BUTTON_RIGHT){
		if (state == GLFW_PRESS){
			cout << "R Click captured by Tool at: " << x << " " << y << endl;
			rClick = true;
			x0 = x;
			y0 = y;
		}
		else{
			rClick = false;
		} 

	}
}


void Tool::onClick(int button, int state, int mods){
	captureClick(button, state, mods);		
}


void Tool::onMouseMove(int x, int y){

}


void Tool::onScroll(double dx, double dy){

}


void onClick1(GLFWwindow* window, int button, int state, int mods){
	if (!Tool::activeTools.empty())	Tool::activeTools.front()->onClick(button, state, mods);
}


void onMouseMove1(GLFWwindow* window, double x, double y){
	if (!Tool::activeTools.empty())	Tool::activeTools.front()->onMouseMove(x, y);
}

void onScroll1(GLFWwindow* window, double dx, double dy){
	if (!Tool::activeTools.empty())	Tool::activeTools.front()->onScroll(dx, dy);
}







	




