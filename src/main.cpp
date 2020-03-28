#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <cmath>
#include <algorithm>
using namespace std;

#include"../include/quickgl.h"


#include "../include/pointcloud.h"


void processInput(GLFWwindow *window);


vector <int> y_slices(float * pos, int n, float res){
	sort_by_y(pos, pos+3*n); 
	vector <int> slice_indices(1,0);
	float min = int(pos[3*0+1]/res)*res;
	for (int i=0; i< n; ++i){
		if (pos[3*i+1] > (min + res*slice_indices.size()) ){
			slice_indices.push_back(i);
			cout << "slice: " << i << ", value = " << pos[3*i+1] <<  endl;
		}
	}
//	if (*slice_indices.end() != n-1) slice_indices.push_back(n-1);
	return slice_indices;
}

vector <int> z_slices(float * pos, int n, float res){
	sort_by_z(pos, pos+3*n); 
	vector <int> slice_indices(1,0);
	float min = int(pos[3*0+2]/res)*res;
	for (int i=0; i< n; ++i){
		if (pos[3*i+2] > (min + res*slice_indices.size()) ){
			slice_indices.push_back(i);
			cout << "slice: " << i << ", value = " << pos[3*i+2] <<  endl;
		}
	}
//	if (*slice_indices.end() != n-1) slice_indices.push_back(n-1);
	return slice_indices;
}


class Particle{
	public:
	float x,y,z;
	int gid;
	int cluster_size;
};

int main(int argc, char **argv){

	cout << "----> Read Data" << endl;
	PointCloud cr;
	cr.read_las("/Users/sachita/Documents/treeseg_project/2017-02-20_21-47-24.las");

	cout << "----> Create DEM" << endl;
	cr.createDEM(0.5,0.5);
//	//cr.dem.printToFile("dem.txt");

	cout << "----> Subtract DEM" << endl;
	cr.subtractDEM_bil();
	cr.deleteGround(0.4); // argument is the height threshold below which to delete points. if understory is thick, it will be difficult to delete ground
//	cr.generateRandomClusters(1000000, 500, -100, 100, -100, 100, 0, 10, 2);

//	cr.write_las("data/2017-02-20_21-47-24_processed.las");
	
	Palette p(1000000);
	//p.createRainbow();
	p.createRandom();

	cout << "----> Group... " << endl;
	cr.group_grid_fof64(0.02);
//	cr.group_grid_hashSTL(0.05);

	cout << "----> Denoise" << endl;	

	// count connectivity
	vector <int> cluster_size(cr.nverts, 0);
	for (int i=0; i<cr.nverts; ++i){
		++cluster_size[cr.group_ids[i]];
	}
	
	// remove extremely small groups
	for (int i=0; i<cr.nverts; ++i){
		if (cluster_size[cr.group_ids[i]] < 3){
			cr.points[3*i+0] = cr.points[3*i+1] = cr.points[3*i+2] = 0;
//			cr.group_ids[i] = cr.nverts;
		}
	}


	// sort by property
	vector <int> sort_idx(cr.nverts);
	for (int i=0; i<cr.nverts; ++i) sort_idx[i]=i;
	sort(sort_idx.begin(), sort_idx.end(), [&cr](int a, int b){return cr.group_ids[a] < cr.group_ids[b];});
	
	
	vector <float> pos2(cr.nverts*3);
	vector <float> gids2(cr.nverts);
	vector <float> clsiz(cr.nverts);
	for (int i=0; i<cr.nverts; ++i){
		pos2[3*i+0] = cr.points[3*sort_idx[i]+0];
		pos2[3*i+1] = cr.points[3*sort_idx[i]+1];
		pos2[3*i+2] = cr.points[3*sort_idx[i]+2];
		gids2[i] = cr.group_ids[sort_idx[i]];
		clsiz[i] = log(1+cluster_size[gids2[i]]);	// log-ed for coloring
	}
	

	cout << "----> Draw" << endl;	
	vector <float> cols9z = p.mapValues(gids2.data(), cr.nverts, 1, 0);	// map group ID
//  vector <float> cols9z = p.mapValues(clsiz.data(), cr.nverts, 1, 0);	// map Cluster size
//	vector <float> cols9z = p.mapValues(cr.points.data(), cr.nverts, 3, 2, 0);	// map z

	

	initQuickGL();


	// Vectors to sample points randomly and store them
	vector <float> vertices_sampled(3*cr.nverts/10);
	vector <float> gids_sampled(cr.nverts/10);

	for(int iter=0; iter<cr.nverts/10; iter++){

		vertices_sampled[3*iter+0] = pos2[3*iter*10+0];
		vertices_sampled[3*iter+1] = pos2[3*iter*10+1];
		vertices_sampled[3*iter+2] = pos2[3*iter*10+2];

		gids_sampled[iter] = gids2[iter*10];

	}

	// Map gid to colour
	vector <float> colours_sampled = p.mapValues(gids_sampled.data(), cr.nverts/10, 1, 0);	// map group ID

	Shape pt(cr.nverts/10, GL_POINTS); //, 4, -1, 1);
	pt.setVertices(vertices_sampled.data());	
	pt.setColors(colours_sampled.data());
	pt.autoExtent();
 
 
//	vector <int> slices = z_slices(cr.points.data(), cr.nverts, 1);
////	vector <int> slices = y_slices(cr.points.data(), cr.nverts, 10);
//	for (int i=0; i<slices.size(); ++i){
//		cout << "slices: " << slices[i] << ": " << cr.points[3*slices[i]+2] << "\n";
//	}
//	cout << endl;

//	vector <Shape*> slices_shapes(slices.size());
//	for (int i=1; i<slices.size(); ++i){
//		int nv = slices[i]-slices[i-1];
//		slices_shapes[i-1] = new Shape(nv, 3, "points", false);
//		slices_shapes[i-1]->pointSize = 1;
//		slices_shapes[i-1]->setVertices(&cr.points[3*slices[i-1]]);
//		slices_shapes[i-1]->setColors(&cols9z[4*slices[i-1]]);
//		slices_shapes[i-1]->setExtent(ex);
//	}
//	slices_shapes[0]->b_render = true;

	
	int current_render = 0;

	float pos3[] = {-20,0,0, 100,0,0, 0,-20,0, 0,100,0, 0,0,-20, 0,0,100};
	float col3[] = {1,0,0, 0.5	,
					 1,0,0, 0.5,
					 0,1,0, 0.5,
					 0,1,0, 0.5,
					 0.0,0.8,1, 0.5,
					 0.0,0.8,1, 0.5
				   };
	
	Shape axis(6, GL_LINES);
	axis.setVertices(pos3);
	axis.setColors((float*)col3);
	

	Camera cam(glm::vec3(2.f,2.f,2.f), glm::vec3(0.5f, 0.5f, 0.f), glm::vec3(0.0f, 0.0f, 1.0f));
	cam.activate();

	// CameraController c;

	
//	while(1){	   // infinite loop needed to poll anim_on signal.
//		slices_shapes[current_render]->b_render = false;
//		current_render = generic_count % (slices_shapes.size()-1);
////		cout << "Currently rendering: " << current_render << "\n";
//		slices_shapes[current_render]->b_render = true;
		
		// glutMainLoop();
//		usleep(20000);
//	}
	// launch sim end.

	while (!glfwWindowShouldClose(window))
    {
        // input
        // -----
        processInput(window);

        // render
        // ------
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

//        // draw our first triangle
//        glUseProgram(shaderProgram);
//        glBindVertexArray(VAO); // seeing as we only have a single VAO there's no need to bind it every time, but we'll do so to keep things a bit more organized
//        //glDrawArrays(GL_TRIANGLES, 0, 6);
//        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0); //points: GL_POINTS
//        // glBindVertexArray(0); // no need to unbind it every time 
 		pt.render();
 		axis.render();
 		
        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // -------------------------------------------------------------------------------
        glfwSwapBuffers(window);
        //glfwPollEvents();
        glfwWaitEvents();
    }


	return 0;
}

// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------
void processInput(GLFWwindow *window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}



