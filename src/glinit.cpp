#include "../include/glinit.h"

#include <iostream>

// #include "../include/shape.h"
// #include "../include/camera.h"
// #include "../include/tool_base.h"

using namespace std;

// GLuint vao;
int width = 800;
int height = 600;

GLFWwindow* window = nullptr;

string errorString(GLenum err){
	if (err == GL_INVALID_ENUM) return "GL_INVALID_ENUM";
	else if (err == GL_INVALID_VALUE) return "GL_INVALID_VALUE";
	else if (err == GL_INVALID_OPERATION) return "GL_INVALID_OPERATION";
	else if (err == GL_STACK_OVERFLOW) return "GL_STACK_OVERFLOW";
	else if (err == GL_STACK_UNDERFLOW) return "GL_STACK_UNDERFLOW";
	else if (err == GL_OUT_OF_MEMORY) return "GL_OUT_OF_MEMORY";
	else return "NO ERROR"; 
//	else if (err == GL_TABLE_TOO_LARGE) return "GL_TABLE_TOO_LARGE";	
}

// // TODO: Should there be a GLInterface class that controls the globals? 
// //       e.g. It will store the active shapes list and active tools list, plus active window, active tool, active camera etc. 
// //            Can also separate input form processing
// // 		      Maybe not?
// //            Maybe yes!

void checkGLError(const char * file, int line){
	GLenum err;
	while((err = glGetError()) != GL_NO_ERROR){
		cout << file << ":" << line << "   Error: " << errorString(err) << endl;
	}
}


void printStatus(const char *step, GLuint context, GLuint status){
	GLint result = GL_FALSE;
	CHECK_GL_ERROR();
	glGetShaderiv(context, status, &result);
	CHECK_GL_ERROR();
	if (result == GL_FALSE) {
		char buffer[1024];
		if (status == GL_COMPILE_STATUS)
			glGetShaderInfoLog(context, 1024, NULL, buffer);
		else
			glGetProgramInfoLog(context, 1024, NULL, buffer);
		if (buffer[0])
			fprintf(stderr, "%s: %s\n", step, buffer);
	}
}


// bool windowIsOpen(){
// 	return !glfwWindowShouldClose(window);
// }

// void onDisplay(void){

//         // render
//         // ------
// //        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
// 		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
// 		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

// 		for (auto it : Shape::allShapes) it->render();	
 
//         // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
//         // -------------------------------------------------------------------------------
//         glfwSwapBuffers(window);

// }


// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void onResize(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and 
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}


// void onClick(GLFWwindow* window, int button, int state, int mods){
// 	if (!Tool::activeTools.empty())	Tool::activeTools.front()->onClick(button, state, mods);
// }


// void onMouseMove(GLFWwindow* window, double x, double y){
// 	if (!Tool::activeTools.empty())	Tool::activeTools.front()->onMouseMove(x, y);
// }

// void onScroll(GLFWwindow* window, double dx, double dy){
// 	if (!Tool::activeTools.empty())	Tool::activeTools.front()->onScroll(dx, dy);
// }

int initQuickGL(){

    // glfw: initialize and configure
    // ------------------------------
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // uncomment this statement to fix compilation on OS X
#endif

    // glfw window creation
    // --------------------
    
    window = glfwCreateWindow(width, height, "LearnOpenGL", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, onResize);

//    // gl3w: load all OpenGL function pointers
//    // ---------------------------------------
	if (gl3wInit()) {
		printf("failed to initialize OpenGL\n");
		return -1;
	}

    // To enable point sizes 
    glEnable(GL_PROGRAM_POINT_SIZE);


	// BUFFERS etc

	// glGenVertexArrays(1, &vao);
	// glBindVertexArray(vao);

	// callbacks
    // glfwSetFramebufferSizeCallback(window, onResize); //CHANGE LATER?
//    glfwSetMouseButtonCallback(window, onClick);
//	glfwSetCursorPosCallback(window, onMouseMove);
//	glfwSetScrollCallback(window, onScroll);
//	glutMouseFunc(onClick);
//	glutMotionFunc(onMouseMove);
	return 0;
}


// void closeQuickGL(){
// 	glBindVertexArray(0);
// 	glDeleteVertexArrays(1, &vao);
// 	glfwTerminate();

// }

