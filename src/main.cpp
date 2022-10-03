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


class Particle{
	public:
	float x,y,z;
	int gid;
	int gsize;
	int parent;
	int rank;
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Union-Find functions to calculate group Indices 
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// given array of parents, find the root of q
int root(int q, Particle * pvec){
	while (q != pvec[q].parent){					   // while (q != par[q])
		pvec[q].parent = pvec[pvec[q].parent].parent;  // par[q] = par[par[q]]  // perform path-compression on the go
		q = pvec[q].parent;							   // q = par[q]            
	}
	return q;
}

// check if p and q have same root, i.e. belong to same group
bool find(int p, int q, Particle * pvec){
	return root(p, pvec) == root(q, pvec);
}

// put p and q in same group by merging the smaller tree into large one
void unite(int p, int q, Particle * pvec){
	int i = root(p, pvec);
	int j = root(q, pvec);
	if (i==j) return;	// if both already have the same root do nothing
	if (pvec[i].rank < pvec[j].rank) {pvec[i].parent = j; pvec[j].rank += pvec[i].rank;}
	else {pvec[j].parent = i; pvec[i].rank += pvec[j].rank; }
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


int main(int argc, char **argv){

	cout << "----> Read Data" << endl;
	PointCloud cr;
	cr.read_las("data/2017-02-20_21-47-24.las");

	cout << "----> Create DEM" << endl;
	cr.createDEM(0.5,0.5);
//	//cr.dem.printToFile("dem.txt");

	cout << "----> Subtract DEM" << endl;
	cr.subtractDEM_bil();
	cr.deleteGround(0.4); // argument is the height threshold below which to delete points. if understory is thick, it will be difficult to delete ground
//	cr.generateRandomClusters(1000000, 500, -100, 100, -100, 100, 0, 10, 2);

//	cr.write_las("data/2017-02-20_21-47-24_processed.las");
	
	Palette p(1000000);
//	p.createRainbow();
	p.createRandom();

	cout << "----> Group... " << endl;
	cr.group_grid_fof64(0.02);
//	cr.group_grid_hashSTL(0.05);


	// -- POSTPROCESS --	
	vector <Particle> pvec(cr.nverts);

	// 1. sort by group id
	vector <int> sort_idx(cr.nverts);
	vector <int> gids2(cr.nverts);
	for (int i=0; i<cr.nverts; ++i) sort_idx[i]=i;
	sort(sort_idx.begin(), sort_idx.end(), [&cr](int a, int b){return cr.group_ids[a] < cr.group_ids[b];});

	// save sorted gid array for easy access for now
	for (int i=0; i<cr.nverts; ++i){
		gids2[i] = cr.group_ids[sort_idx[i]];
	}

	// 2. initialize ranks and parents and sorted points
	for (int i=0; i<cr.nverts; ++i){
		pvec[i].x = cr.points[3*sort_idx[i]+0];
		pvec[i].y = cr.points[3*sort_idx[i]+1];
		pvec[i].z = cr.points[3*sort_idx[i]+2];
		pvec[i].parent = i;
		pvec[i].rank = 1;
	}
	
	// 3. re-construct grouping tree
	for (int j=0; j<pvec.size()-1; ++j){
		if (gids2[j+1] == gids2[j]){
			unite(j,j+1, pvec.data());
		}
	}
	
	// 4. reassign groups
	for (int j=0; j<pvec.size(); ++j){
		pvec[j].gid = root(j, pvec.data());
	}

	// 5. count connectivity
	vector <int> cluster_size(pvec.size(), 0);
	for (int i=0; i<pvec.size(); ++i){
		++cluster_size[pvec[i].gid];
	}
	for (int i=0; i<pvec.size(); ++i){
		pvec[i].gsize = cluster_size[pvec[i].gid];
	}

	
	// output processed cloud to file
	ofstream fout("data/processed_particles.txt");
	fout <<  "x" << "\t" << "y" << "\t" << "z" << "\t"
		 <<  "gid" << "\t"
		 <<  "gsize" << "\t"
		 <<  "parent" << "\t"
		 <<  "rank" << "\n";
	
	for (auto p : pvec){
		fout <<  p.x << "\t" << p.y << "\t" << p.z  << "\t"
			 <<  p.gid << "\t"
			 <<  p.gsize << "\t"
			 <<  p.parent << "\t"
			 <<  p.rank << "\n";
	}
	fout.close();
		
	// -----------------
	


	cout << "----> Denoise" << endl;	
	
	// remove extremely small groups
	for (int i=0; i<pvec.size(); ++i){
		if (pvec[i].gsize < 20){
			pvec[i].x = pvec[i].y = pvec[i].z = 0;
		}
	}


	cout << "----> Draw" << endl;	

	vector<float> gids3(pvec.size());
	vector<float> pos2(3*pvec.size());
	for (int i=0; i<pvec.size(); ++i){
		gids3[i] = pvec[i].gid;
		pos2[3*i+0] = pvec[i].x;
		pos2[3*i+1] = pvec[i].y;
		pos2[3*i+2] = pvec[i].z;
	}
	
	vector <float> cols9z = p.mapValues(gids3.data(), cr.nverts, 1, 0);	// map group ID
//	vector <float> cols9z = p.mapValues(clsiz.data(), cr.nverts, 1, 0);	// map Cluster size
//	vector <float> cols9z = p.mapValues(cr.points.data(), cr.nverts, 3, 2, 0);	// map z


	initQuickGL(argc, argv);


	Shape pt(cr.nverts, GL_POINTS); //, 4, -1, 1);
	pt.setVertices(pos2.data());	
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


