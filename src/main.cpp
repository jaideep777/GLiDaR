#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <cmath>
#include <algorithm>
using namespace std;

#include <quickgl.h>


#include "../headers/pointcloud.h"


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




int main(int argc, char **argv){

	initQuickGL(argc, argv);

	PointCloud cr;
	cr.read_las("data/2017-02-20_21-47-24.las");
	cr.createDEM(0.5,0.5);
//	//cr.dem.printToFile("dem.txt");
	cr.subtractDEM_bil();
	cr.deleteGround(0.4); // argument is the height threshold below which to delete points. if understory is thick, it will be difficult to delete ground
//	cr.generateRandomClusters(1000000, 500, -100, 100, -100, 100, 0, 10, 2);
	cr.write_ascii("data/2017-02-20_21-47-24.txt");
	
	Palette p(1000000);
	p.createRainbow();
//	p.createRandom();


//	vector <float> cols9z = p.mapValues(gids.data(), cr.nverts, 1, 0);	// map group ID
	vector <float> cols9z = p.mapValues(cr.points.data(), cr.nverts, 3, 2, 0);	// map z

	Shape pt(cr.nverts, GL_POINTS); //, 4, -1, 1);
	pt.setVertices(cr.points.data());	
	pt.setColors(&cols9z[0]);
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

	CameraController c;

	
//	while(1){	   // infinite loop needed to poll anim_on signal.
//		slices_shapes[current_render]->b_render = false;
//		current_render = generic_count % (slices_shapes.size()-1);
////		cout << "Currently rendering: " << current_render << "\n";
//		slices_shapes[current_render]->b_render = true;
		
		glutMainLoop();
//		usleep(20000);
//	}
	// launch sim end.
	
	return 0;
}


