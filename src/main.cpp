#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <cmath>
#include <algorithm>
using namespace std;

#include <quickgl.h>

#include "../utils/simple_math.h"
#include "../utils/simple_timer.h"

#include "../headers/hashtable.h"
#include "../headers/pointcloud.h"

#include <unordered_map>



Grid par;


vector <int> y_slices(float * pos, int n, float res){
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

	PointCloud cr;
//	cr.read_las("/home/chethana/Data/LIDAR/2017-02-20_21-47-24.las");
//	cr.createDEM(0.5,0.5);
//	//cr.dem.printToFile("dem.txt");
//	cr.subtractDEM();
////	cr.deleteGround(0.2);
	cr.generateRandomClusters(1000000, 500, -100, 100, -100, 100, 0, 10, 2);
	
	initQuickGL(argc, argv);

	Palette p(1000000);
//	p.createRainbow();
	p.createRandom();
	
	cout << "sort...\n";
	sort_by_z(cr.points.data(), cr.points.data()+3*cr.nverts); 
	for (int i=0; i<10; ++i) cout << cr.points[3*i] << " " << cr.points[3*i+1] << " " << cr.points[3*i+2] << endl;

//	cr.group_serial(2);
//	cr.denoise(2);
//	cr.group_grid_hash(0.2);
	cr.group_grid_hash2(2);
//	cr.group_grid_hashSTL(2);
//	cr.group_grid_map(0.5);
//	cr.denoise(1);
//	cr.countNeighbours_hash_gpu(1);


	Camera cam(glm::vec3(2.f,2.f,2.f), glm::vec3(0.5f, 0.5f, 0.f), glm::vec3(0.0f, 0.0f, 1.0f));
	cam.activate();


	CameraController c;

//	vector <float> cols9z = p.map_values(&cr.points[2], cr.nverts, 3);	// map z value
	vector <float> gids(cr.group_ids.begin(), cr.group_ids.end());
	vector <float> cols9z = p.mapValues(gids.data(), cr.nverts, 1, 0);	// map group ID
//	vector <float> nn(cr.neighbourCounts.begin(), cr.neighbourCounts.end());
//	vector <float> cols9z = p.map_values(nn.data(), cr.nverts, 1);	// map neighbour counts
	Shape pt(cr.nverts, GL_POINTS); //, 4, -1, 1);
	pt.setVertices(cr.points.data());	
	pt.setColors(&cols9z[0]);
//	vector <float> ex = calcExtent(cr.points.data(), cr.nverts, 3);
	pt.autoExtent();
//	pt.pointSize = 2;
 
 
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




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Union-Find functions to calculate group Indices 
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// given array of parents, find the root of q
int root(int q, int* par){
	while (q != par[q]){
		par[q] = par[par[q]];
		q = par[q];
	}
	return q;
}

// check if p and q have same root, i.e. belong to same group
bool find(int p, int q, int *par){
	return root(p, par) == root(q, par);
}

// put p and q in same group by merging the smaller tree into large one
void unite(int p, int q, int *par, int *sz){
	int i = root(p, par);
	int j = root(q, par);
	if (i==j) return;	// if both already have the same root do nothing
	if (sz[i] < sz[j]) {par[i]=j; sz[j] += sz[i];}
	else 			   {par[j]=i; sz[i] += sz[j];}
}


void PointCloud::group_serial(float Rg){
	
	SimpleTimer T; T.start();
	
	vector <int> par(nverts);
	vector <int> sz(nverts,1);
	for (int i=0; i<nverts; ++i) par[i] = i;
	
	long long int pairs = 0;
	
	for (int p=0; p<nverts; ++p){
		for (int q=0; q<= p; ++q) {

			float d2other = distance(p, q);
//			if (p % 100 == 0)
//			cout << p << "," << q << ": " << d2other << endl;
			
			// if distance is < R_grp, assign same group
			if (d2other < Rg){
				unite(p,q, par.data(), sz.data());
			} 
			++pairs;
		}
	}

	group_ids.resize(nverts);
	for (int i=0; i<nverts; ++i){
		group_ids[i] = root(i, par.data());
	}

	T.stop(); T.printTime();
	
	cout << "Pairs compared = " << pairs << endl;
//	for (int i=0; i<group_ids.size(); ++i){
//		cout << sz[i] << " ";
//	}
//	cout << endl;
	
}








int3 cellIndex(float3 pos, float3 origin){
	int3 cell;
	cell.x = (pos.x - origin.x)/par.cellSize.x;
	cell.y = (pos.y - origin.y)/par.cellSize.y;
	cell.z = (pos.z - origin.z)/par.cellSize.z;
	return cell;
}

unsigned int cellHash(int3 cell){
	//		   iz * nx * ny                     +    iy * nx            + ix
	return cell.z*par.gridSize.x*par.gridSize.y + cell.y*par.gridSize.x + cell.x;
	// can use z-order curve as hash here rather than 1D cell-index
}

int clamp(int x, int xmin, int xmax){
	if (x < xmin) x = xmin;
	if (x > xmax) x = xmax;
	return x;
}


void PointCloud::group_grid(float Rg){
	
	SimpleTimer T; T.reset(); T.start();
	calcGridParams(Rg, par);
	T.stop(); T.printTime("calcGrid");

	
	float3 * pos = (float3*)points.data();

	T.reset(); T.start();	
	vector <unsigned int> point_hashes(nverts);
	vector <unsigned int> point_ids(nverts);
	// get the cell ID for each particle
	for (int i=0; i<nverts; ++i) {
		int3 cell_id = cellIndex(pos[i], par.origin);
		point_hashes[i] = cellHash(cell_id);
		point_ids[i]    = i;
	}
	T.stop(); T.printTime("cell Ids");

//	cout << "points, hashes: " << endl;	
//	for (int i=0; i<10; ++i) cout << point_ids[i] << " " << point_hashes[i] << "\n";

	// sort particles by cell ID
	T.reset(); T.start();	
	sort(point_ids.begin(), point_ids.end(), [&point_hashes](unsigned int i, unsigned int j){return point_hashes[i] < point_hashes[j];}); 
	T.stop(); T.printTime("sort");
	
//	cout << "points, hashes: " << endl;	
//	for (int i=0; i<10; ++i) cout << point_ids[i] << " " << point_hashes[point_ids[i]] << "\n";

	T.reset(); T.start();
	// calc start and end of each cell (both reflect indices in sorted hashes list*)
	//  -- *sorted hashes list is accessed as point_hashes[point_ids[0:nverts]]
	int ncells = par.gridSize.x*par.gridSize.y*par.gridSize.z;
	vector <int> cell_start(ncells, -1);
	vector <int> cell_end  (ncells, -1);
	for (int i=1; i<nverts; ++i){
		int k = point_ids[i];		// get ith point in sorted list
		int kprev = point_ids[i-1];	// i-1 th point in sorted list
		
		if (i==1) {
			cell_start[point_hashes[kprev]] = i-1;
		}
		if (i==nverts-1){
			cell_end[point_hashes[k]] = i;
		}

		if (point_hashes[k] != point_hashes[kprev]){
			cell_start[point_hashes[k]] = i;
			cell_end[point_hashes[kprev]] = i-1;
		}
	}
	T.stop(); T.printTime("cell s/e");
//	cout << "ALL:" << endl;
//	for (int i=0; i<nverts; ++i)
//	cout << point_ids[i] << " " << point_hashes[point_ids[i]] << " " << endl;
//	
//	cout << "CELLS" << endl;
//	for (int i=0; i<ncells; ++i)
//	cout << i << " " << cell_start[i] << " " << cell_end[i] << endl;
	
	vector <int> parents(nverts);
	vector <int> sz(nverts,1);
	for (int i=0; i<nverts; ++i) parents[i] = i;
	
	long long int pairs = 0;
	
	T.reset(); T.start();
	// For each particle, check particles in all neighbouring grids and group them
	for (int i=0; i<nverts; ++i){
		
		float3 p1 = pos[i];
		int3 cell = cellIndex(p1, par.origin);	// get particle cell
//		cout << "cell: " << cell.x << ", " << cell.y << ", " << cell.z << endl;
		for (int ix = -1; ix <=1; ++ix){
			for (int iy = -1; iy <=1; ++iy){
				for (int iz = -1; iz <=1; ++iz){
					int3 cell_new;
					cell_new.x = clamp(cell.x + ix, 0, par.gridSize.x-1);
					cell_new.y = clamp(cell.y + iy, 0, par.gridSize.y-1);
//					cout << "\tcell: " << cell_new.x << ", " << cell_new.y << ", " << cell_new.z << endl;
					
					unsigned int hash = cellHash(cell_new);
					
					// loop over all particles in the cell
					for (int k=cell_start[hash]; k<=cell_end[hash]; ++k){
						float3 p2 = pos[point_ids[k]];
						
						if (distance(i, point_ids[k]) < Rg) unite(i, point_ids[k], parents.data(), sz.data());
						++pairs;
						
//						if (pairs % 1000000 == 0) cout << "pairs : " << pairs << endl;
					}
				}
			}
		}
		
	}
	T.stop(); T.printTime("pairs");
	
	T.reset(); T.start();
	group_ids.resize(nverts);
	for (int i=0; i<nverts; ++i){
		group_ids[i] = root(i, parents.data());
	}
	T.stop(); T.printTime("gid");

	
	cout << "Pairs compared = " << pairs << endl;

	
}


int particles_compared_per_cellpair = 0;
int count_pcpc = 0;

int3 PointCloud::mergeCells(int c1, int c2, float Rg, map <int, int2> & cells, vector <unsigned int> & pt_ids){
	
	int count = 0;
	
	int3 result; // (yes/no, point_id1, point_id2)
	result.x = result.y = result.z = 0;

	map <int, int2>::iterator ic1 = cells.find(c1);
	map <int, int2>::iterator ic2 = cells.find(c2);
	if (ic1 == cells.end() || ic2 == cells.end()) return result; // if either cell doesnt exist, cant merge cells
	
	for (int p1 = ic1->second.x; p1 <= ic1->second.y; ++p1){	
		for (int p2 = ic2->second.x; p2 <= ic2->second.y; ++p2){
			++count;
			
			int ip1 = pt_ids[p1];
			int ip2 = pt_ids[p2];
			if (distance(ip1, ip2) < Rg){
				result.x = 1;
				result.y = ip1;
				result.z = ip2;
				
				particles_compared_per_cellpair += count;
				count_pcpc += 1;
				return result;
			} 
		}
	}
	
	return result;
}


int cells_compared_per_cell = 0;
int count_ccpc = 0;

void PointCloud::group_grid_map(float Rg){

	SimpleTimer T; T.reset(); T.start();
	calcGridParams(Rg/sqrt(3), par);
	T.stop(); T.printTime("calcGrid");

	
	float3 * pos = (float3*)points.data();

	T.reset(); T.start();	
	vector <unsigned int> point_hashes(nverts);
	vector <unsigned int> point_ids(nverts);
	// get the cell ID for each particle
	for (int i=0; i<nverts; ++i) {
		int3 cell_id = cellIndex(pos[i], par.origin);
		point_hashes[i] = cellHash(cell_id);
		point_ids[i]    = i;
	}
	T.stop(); T.printTime("cell Ids");


	// sort particles by cell ID
	T.reset(); T.start();	
	sort(point_ids.begin(), point_ids.end(), [&point_hashes](unsigned int i, unsigned int j){return point_hashes[i] < point_hashes[j];}); 
	T.stop(); T.printTime("sort");

	vector <int> parents(nverts);
	vector <int> sz(nverts,1);
	for (int i=0; i<nverts; ++i) parents[i] = i;

//	cout << "ALL:" << endl;
//	for (int i=0; i<30; ++i)
//	cout << point_ids[i] << " " << point_hashes[point_ids[i]] << " " << endl;

	map <int, int2> filled_cells;
	
	// Merge all points in a given cell into same group
	T.reset(); T.start();
	// start and end of each cell (both reflect indices in sorted hashes list*)
	//  -- *sorted hashes list is accessed as point_hashes[point_ids[0:nverts]]
	int start = 0, next=1;
	int cellStart = point_ids[start];
	while (next < nverts){
		int cellNext = point_ids[next];
		if (point_hashes[cellNext] == point_hashes[cellStart]){
			unite(cellNext, cellStart, parents.data(), sz.data());
			//cout << "uniting " << cellStart << "|" << cellNext << " (" << point_hashes[cellStart] << "|" << point_hashes[cellNext] << endl;  
		}
		else {
			int2 se; se.x = start; se.y = next-1;
			filled_cells[point_hashes[cellStart]] = se;

			start = next;
			cellStart = cellNext; 
		}
		++next;
	}
	T.stop(); T.printTime("cell s/e");


//	cout << "ALL:" << endl;
//	for (int i=0; i<30; ++i)
//	cout << point_ids[i] << " " << point_hashes[point_ids[i]] << " " << group_ids[point_ids[i]] << endl;

//	cout << "MAP (" << filled_cells.size() << "):" << endl;	
//	int i=0;
//	for (map<int,int2>::iterator it = filled_cells.begin(); it != filled_cells.end(); ++it){
//		cout << it->first << " " << it->second.x << " " << it->second.y << "\n";
//		++i; if (i >10) break;
//	}
//	cout << (--filled_cells.end())->first << " " << (--filled_cells.end())->second.x << " " << (--filled_cells.end())->second.y << endl;
//	// 
	
	cout << "No. of filled cells = " << filled_cells.size() << endl;
	
	// group grid cells
	long long int pairs = 0;
	T.reset(); T.start();
	for (map<int,int2>::iterator it = filled_cells.begin(); it != filled_cells.end(); ++it){

		int c1 = it->first;

		int3 cell;
		cell.z = int(c1/par.gridSize.x/par.gridSize.y);
		cell.y = int((c1 - cell.z*par.gridSize.x*par.gridSize.y)/par.gridSize.x);
		cell.x = int((c1 - cell.z*par.gridSize.x*par.gridSize.y - cell.y*par.gridSize.x));

		int count=0;
		
		for (int ix = -2; ix <=2; ++ix){
			for (int iy = -2; iy <=2; ++iy){
				for (int iz = -2; iz <=2; ++iz){
					
					int3 cell_new;
					cell_new.x = clamp(cell.x + ix, 0, par.gridSize.x-1);
					cell_new.y = clamp(cell.y + iy, 0, par.gridSize.y-1);
					cell_new.z = clamp(cell.z + iz, 0, par.gridSize.z-1);
//					cout << "\tcell: " << cell_new.x << ", " << cell_new.y << ", " << cell_new.z << endl;

					int c2 = cellHash(cell_new);
					
					if (c2 < c1){	
						int3 res = mergeCells(c1, c2, Rg, filled_cells, point_ids);
						if (res.x == 1) unite(res.y, res.z, parents.data(), sz.data());

						++count;
						++pairs;
					}
//					if (pairs % 10000 == 0) cout << pairs << "pairs compared.\n"; 
				}
			}
		}
		
		cells_compared_per_cell += count;
		count_ccpc += 1;
	}
	cout << "pairs = " << pairs << endl;
	T.stop(); T.printTime("group cells");
	
	T.reset(); T.start();
	group_ids.resize(nverts);
	for (int i=0; i<nverts; ++i){
		group_ids[i] = root(i, parents.data());
	}
	T.stop(); T.printTime("gid");
	
	cout << "Summary: \n";
	cout << "Particle pairs compared per mergeCell: " <<  float(particles_compared_per_cellpair)/count_pcpc << endl;
	cout << "Cells compared per cell: " <<  float(cells_compared_per_cell)/count_ccpc << endl;


}



#define printCell(cell) "(" << cell.x << "," << cell.y << "," << cell.z << ")"

int3 PointCloud::mergeCells_hash(int3 c1, int3 c2, float Rg, HashNode * ht, vector <unsigned int> & pt_ids, int length){
	
	int count = 0;
	
	int3 result; // (yes/no, point_id1, point_id2)
	result.x = result.y = result.z = 0;
	
	int attempts1, attempts2;
	int id_c1 = hash_find(c1, ht, length, &attempts1);
	int id_c2 = hash_find(c2, ht, length, &attempts2);
	if (id_c1 == -1 || id_c2 == -1) return result; // if either cell doesnt exist, cant merge cells
	
	for (int p1 = ht[id_c1].value.x; p1 <= ht[id_c1].value.y; ++p1){	
		for (int p2 = ht[id_c2].value.x; p2 <= ht[id_c2].value.y; ++p2){
			++count;
			
			int ip1 = pt_ids[p1];
			int ip2 = pt_ids[p2];
			if (distance(ip1, ip2) < Rg){
				result.x = 1;
				result.y = ip1;
				result.z = ip2;
				
				particles_compared_per_cellpair += count;
				count_pcpc += 1;
				return result;
			} 
		}
	}
//	cout << "attempts1 = " << attempts1 << " " << "attempts2 = " << attempts2 << endl;
	return result;
}

void PointCloud::group_grid_hash(float Rg){

	SimpleTimer T; T.reset(); T.start();
	calcGridParams(Rg/sqrt(3), par);
	T.stop(); T.printTime("calcGrid");

	
	float3 * pos = (float3*)points.data();

	T.reset(); T.start();	
	vector <unsigned int> point_hashes(nverts);
	vector <unsigned int> point_ids(nverts);
	// get the cell ID for each particle
	for (int i=0; i<nverts; ++i) {
		int3 cell_id = cellIndex(pos[i], par.origin);
		point_hashes[i] = cellHash(cell_id);
		point_ids[i]    = i;
	}
	T.stop(); T.printTime("cell Ids");


	// sort particles by cell ID
	T.reset(); T.start();	
	sort(point_ids.begin(), point_ids.end(), [&point_hashes](unsigned int i, unsigned int j){return point_hashes[i] < point_hashes[j];}); 
	T.stop(); T.printTime("sort");

	vector <int> parents(nverts);
	vector <int> sz(nverts,1);
	for (int i=0; i<nverts; ++i) parents[i] = i;

//	cout << "ALL:" << endl;
//	for (int i=0; i<30; ++i)
//	cout << point_ids[i] << " " << point_hashes[point_ids[i]] << " " << endl;


	int hashTable_size = 6000011;
	vector <HashNode> hashTable(hashTable_size);

	int n_attempts = 0;	
	int avg_attempts = 0;
	int max_attempts = 0;

	// Merge all points in a given cell into same group
	T.reset(); T.start();
	// start and end of each cell (both reflect indices in sorted hashes list*)
	//  -- *sorted hashes list is accessed as point_hashes[point_ids[0:nverts]]
	int start = 0, next=1;
	int cellStart = point_ids[start];
	while (next < nverts){
		int cellNext = point_ids[next];
		if (point_hashes[cellNext] != point_hashes[cellStart]){
			int2 se; se.x = start; se.y = next-1;
			
			int3 cell;
			cell.z = int(point_hashes[cellStart]/par.gridSize.x/par.gridSize.y);
			cell.y = int((point_hashes[cellStart] - cell.z*par.gridSize.x*par.gridSize.y)/par.gridSize.x);
			cell.x = int((point_hashes[cellStart] - cell.z*par.gridSize.x*par.gridSize.y - cell.y*par.gridSize.x));

			int a = hash_insert(cell, se, hashTable.data(), hashTable_size);
			++n_attempts;
			avg_attempts += a;
			max_attempts = max(max_attempts, a);

			start = next;
			cellStart = cellNext; 
		}
		else {
			unite(cellNext, cellStart, parents.data(), sz.data());
			//cout << "uniting " << cellStart << "|" << cellNext << " (" << point_hashes[cellStart] << "|" << point_hashes[cellNext] << endl;  
		}
		++next;
	}
	cout << "Insertion Attempts: " << float(avg_attempts)/n_attempts << " " << max_attempts << endl;
	T.stop(); T.printTime("map cells");


//	cout << "ALL:" << endl;
//	for (int i=0; i<30; ++i)
//	cout << point_ids[i] << " " << point_hashes[point_ids[i]] << " " << group_ids[point_ids[i]] << endl;

//	cout << "MAP (" << filled_cells.size() << "):" << endl;	
//	int i=0;
//	for (map<int,int2>::iterator it = filled_cells.begin(); it != filled_cells.end(); ++it){
//		cout << it->first << " " << it->second.x << " " << it->second.y << "\n";
//		++i; if (i >10) break;
//	}
//	cout << (--filled_cells.end())->first << " " << (--filled_cells.end())->second.x << " " << (--filled_cells.end())->second.y << endl;
//	// 

	// group grid cells
	int n_filled_cells = 0;
	long long int pairs = 0;
	T.reset(); T.start();
	for (int i=0; i<hashTable_size; ++i){
		
		int3 cell = hashTable[i].key;
		if (cell == HashNode().key) continue;
		++n_filled_cells;
		
		int count=0;
		
		for (int ix = -2; ix <=2; ++ix){
			for (int iy = -2; iy <=2; ++iy){
				for (int iz = -2; iz <=2; ++iz){
					
					int3 cell_new;
					cell_new.x = clamp(cell.x + ix, 0, par.gridSize.x-1);
					cell_new.y = clamp(cell.y + iy, 0, par.gridSize.y-1);
					cell_new.z = clamp(cell.z + iz, 0, par.gridSize.z-1);
//					cout << "\tcell: " << cell_new.x << ", " << cell_new.y << ", " << cell_new.z << endl;

					if (cellHash(cell_new) < cellHash(cell)){	
//						int3 res = mergeCells(c1, c2, Rg, filled_cells, point_ids);
						int3 res = mergeCells_hash(cell, cell_new, Rg, hashTable.data(), point_ids, hashTable_size);
						if (res.x == 1) unite(res.y, res.z, parents.data(), sz.data());

						++count;
						++pairs;
					}
//					if (pairs % 10000 == 0) cout << pairs << "pairs compared.\n"; 
				}
			}
		}
		
		cells_compared_per_cell += count;
		count_ccpc += 1;
	}
	cout << "pairs = " << pairs << endl;
	T.stop(); T.printTime("group cells");
	
	T.reset(); T.start();
	group_ids.resize(nverts);
	for (int i=0; i<nverts; ++i){
		group_ids[i] = root(i, parents.data());
	}
	T.stop(); T.printTime("gid");
	
	cout << "Summary: \n";
	cout << "Particle pairs compared per mergeCell: " <<  float(particles_compared_per_cellpair)/count_pcpc << endl;
	cout << "Cells compared per cell: " <<  float(cells_compared_per_cell)/count_ccpc << endl;
	cout << "No. of filled cells = " << n_filled_cells << endl;
	


}


//class Key
//{	
//	public:
//  int x;
//  int y;
//  int z;

//  bool operator==(const Key &other) const
//  { 
//  	return (x == other.x && y== other.y && z == other.z);
//  }
//};




int3 PointCloud::mergeCells_hash_stl(int3 c1, int3 c2, float Rg, unordered_map<int3,int2,KeyHasher> &ht, vector <unsigned int> & pt_ids, int length){
	
	int count = 0;
	
	int3 result; // (yes/no, point_id1, point_id2)
	result.x = result.y = result.z = 0;
	
	int attempts1, attempts2;
	unordered_map<int3, int2, KeyHasher>::iterator id_c1 = ht.find(c1);
	unordered_map<int3, int2, KeyHasher>::iterator id_c2 = ht.find(c2);
	if (id_c1 == ht.end() || id_c2 == ht.end()) return result; // if either cell doesnt exist, cant merge cells
	
	for (int p1 = id_c1->second.x; p1 <= id_c1->second.y; ++p1){	
		for (int p2 = id_c2->second.x; p2 <= id_c2->second.y; ++p2){
			++count;
			
			int ip1 = pt_ids[p1];
			int ip2 = pt_ids[p2];
			if (distance(ip1, ip2) < Rg){
				result.x = 1;
				result.y = ip1;
				result.z = ip2;
				
				particles_compared_per_cellpair += count;
				count_pcpc += 1;
				return result;
			} 
		}
	}
//	cout << "attempts1 = " << attempts1 << " " << "attempts2 = " << attempts2 << endl;
	return result;
}

void PointCloud::group_grid_hashSTL(float Rg){

	SimpleTimer T; T.reset(); T.start();
	calcGridParams(Rg/sqrt(3), par);
	T.stop(); T.printTime("calcGrid");

	
	float3 * pos = (float3*)points.data();

	T.reset(); T.start();	
	vector <unsigned int> point_hashes(nverts);
	vector <unsigned int> point_ids(nverts);
	// get the cell ID for each particle
	for (int i=0; i<nverts; ++i) {
		int3 cell_id = cellIndex(pos[i], par.origin);
		point_hashes[i] = cellHash(cell_id);
		point_ids[i]    = i;
	}
	T.stop(); T.printTime("cell Ids");


	// sort particles by cell ID
	T.reset(); T.start();	
	sort(point_ids.begin(), point_ids.end(), [&point_hashes](unsigned int i, unsigned int j){return point_hashes[i] < point_hashes[j];}); 
	T.stop(); T.printTime("sort");

	vector <int> parents(nverts);
	vector <int> sz(nverts,1);
	for (int i=0; i<nverts; ++i) parents[i] = i;

//	cout << "ALL:" << endl;
//	for (int i=0; i<30; ++i)
//	cout << point_ids[i] << " " << point_hashes[point_ids[i]] << " " << endl;

//	map <int, int2> filled_cells;
//	
//	// Merge all points in a given cell into same group
//	T.reset(); T.start();
//	// start and end of each cell (both reflect indices in sorted hashes list*)
//	//  -- *sorted hashes list is accessed as point_hashes[point_ids[0:nverts]]
//	int start = 0, next=1;
//	int cellStart = point_ids[start];
//	while (next < nverts){
//		int cellNext = point_ids[next];
//		if (point_hashes[cellNext] != point_hashes[cellStart]){
//			int2 se; se.x = start; se.y = next-1;
//			filled_cells[point_hashes[cellStart]] = se;

//			start = next;
//			cellStart = cellNext; 
//		}
//		else {
//			unite(cellNext, cellStart, parents.data(), sz.data());
//			//cout << "uniting " << cellStart << "|" << cellNext << " (" << point_hashes[cellStart] << "|" << point_hashes[cellNext] << endl;  
//		}
//		++next;
//	}
//	T.stop(); T.printTime("cell s/e");



	unordered_map <int3, int2, KeyHasher> filled_cells_hash(6000011);
	
	// Merge all points in a given cell into same group
	T.reset(); T.start();
	// start and end of each cell (both reflect indices in sorted hashes list*)
	//  -- *sorted hashes list is accessed as point_hashes[point_ids[0:nverts]]
	int start = 0, next=1;
	int cellStart = point_ids[start];
	while (next < nverts){
		int cellNext = point_ids[next];
		if (point_hashes[cellNext] != point_hashes[cellStart]){
			int2 se; se.x = start; se.y = next-1;
			
			int3 cell;
			cell.z = int(point_hashes[cellStart]/par.gridSize.x/par.gridSize.y);
			cell.y = int((point_hashes[cellStart] - cell.z*par.gridSize.x*par.gridSize.y)/par.gridSize.x);
			cell.x = int((point_hashes[cellStart] - cell.z*par.gridSize.x*par.gridSize.y - cell.y*par.gridSize.x));

			pair<int3, int2> pair(cell, se);
			filled_cells_hash.insert(pair);

			start = next;
			cellStart = cellNext; 
		}
		else {
			unite(cellNext, cellStart, parents.data(), sz.data());
			//cout << "uniting " << cellStart << "|" << cellNext << " (" << point_hashes[cellStart] << "|" << point_hashes[cellNext] << endl;  
		}
		++next;
	}
	T.stop(); T.printTime("cell s/e");

//	cout << "Test MAP (" << filled_cells.size() << ")" << endl;
//	cout << "Test UNORDENRE MAP (" << filled_cells_hash.size() << ")" << endl;


//	for (map<int,int2>::iterator it = filled_cells.begin(); it != filled_cells.end(); ++it){
//		int3 cell;
//		cell.z = int(it->first/par.gridSize.x/par.gridSize.y);
//		cell.y = int((it->first - cell.z*par.gridSize.x*par.gridSize.y)/par.gridSize.x);
//		cell.x = int((it->first - cell.z*par.gridSize.x*par.gridSize.y - cell.y*par.gridSize.x));
//		
//		if (filled_cells_hash.find(cell) == filled_cells_hash.end()) cout << "not found\n";
//		else cout << printCell(cell) << " | " << printCell(filled_cells_hash.find(cell)->first) << endl;
//	}

//	for (unordered_map<int3,int2,KeyHasher>::iterator it = filled_cells_hash.begin(); it != filled_cells_hash.end(); ++it){
//		cout << printCell(it->first) << endl;
//	}
	
//	int3 test; test.x = 58; test.y= 92; test.z = 6;
//	if (filled_cells_hash.find(test) == filled_cells_hash.end()) cout << "not found\n";
//	else cout << printCell(test) << " | " << printCell(filled_cells_hash.find(test)->first) << endl;
	

//	cout << "ALL:" << endl;
//	for (int i=0; i<30; ++i)
//	cout << point_ids[i] << " " << point_hashes[point_ids[i]] << " " << group_ids[point_ids[i]] << endl;

//	cout << "MAP (" << filled_cells.size() << "):" << endl;	
//	int i=0;
//	for (map<int,int2>::iterator it = filled_cells.begin(); it != filled_cells.end(); ++it){
//		cout << it->first << " " << it->second.x << " " << it->second.y << "\n";
//		++i; if (i >10) break;
//	}
//	cout << (--filled_cells.end())->first << " " << (--filled_cells.end())->second.x << " " << (--filled_cells.end())->second.y << endl;
//	// 

//	int hashTable_size = 6000011;
//	vector <HashNode> hashTable(hashTable_size);
//	
//	int avg_attempts = 0;
//	int max_attempts = 0;
//	T.reset(); T.start();
//	for (map <int,int2>::iterator it = filled_cells.begin(); it != filled_cells.end(); ++it){
//		
//		int c1 = it->first;
//		int3 cell;
//		cell.z = int(c1/par.gridSize.x/par.gridSize.y);
//		cell.y = int((c1 - cell.z*par.gridSize.x*par.gridSize.y)/par.gridSize.x);
//		cell.x = int((c1 - cell.z*par.gridSize.x*par.gridSize.y - cell.y*par.gridSize.x));

//		int2 val; val.x = it->second.x; val.y = it->second.y;
//		int a = hash_insert(cell, val, hashTable.data(), hashTable_size);
//		avg_attempts += a;
//		max_attempts = max(max_attempts, a);
//		
////		cells.push_back(cell);
//	}
//	cout << "Insertion Attempts: " << float(avg_attempts)/filled_cells.size() << " " << max_attempts << endl;
//	T.stop(); T.printTime("map cells");
//	
//	
	// group grid cells
	long long int pairs = 0;
	T.reset(); T.start();
	for (unordered_map<int3,int2,KeyHasher>::iterator it = filled_cells_hash.begin(); it != filled_cells_hash.end(); ++it){

		int3 c1 = it->first;
		int3 cell = c1;
//		int3 cell;
//		cell.z = int(c1/par.gridSize.x/par.gridSize.y);
//		cell.y = int((c1 - cell.z*par.gridSize.x*par.gridSize.y)/par.gridSize.x);
//		cell.x = int((c1 - cell.z*par.gridSize.x*par.gridSize.y - cell.y*par.gridSize.x));

		int count=0;
		
		for (int ix = -2; ix <=2; ++ix){
			for (int iy = -2; iy <=2; ++iy){
				for (int iz = -2; iz <=2; ++iz){
					
					int3 cell_new;
					cell_new.x = clamp(cell.x + ix, 0, par.gridSize.x-1);
					cell_new.y = clamp(cell.y + iy, 0, par.gridSize.y-1);
					cell_new.z = clamp(cell.z + iz, 0, par.gridSize.z-1);
//					cout << "\tcell: " << cell_new.x << ", " << cell_new.y << ", " << cell_new.z << endl;

					int c2 = cellHash(cell_new);
					
					if (c2 < cellHash(c1)){	
//						int3 res = mergeCells(c1, c2, Rg, filled_cells, point_ids);
						int3 res = mergeCells_hash_stl(cell, cell_new, Rg, filled_cells_hash, point_ids, 0);
						if (res.x == 1) unite(res.y, res.z, parents.data(), sz.data());

						++count;
						++pairs;
					}
//					if (pairs % 10000 == 0) cout << pairs << "pairs compared.\n"; 
				}
			}
		}
		
		cells_compared_per_cell += count;
		count_ccpc += 1;
	}
	cout << "pairs = " << pairs << endl;
	T.stop(); T.printTime("group cells");
	
	T.reset(); T.start();
	group_ids.resize(nverts);
	for (int i=0; i<nverts; ++i){
		group_ids[i] = root(i, parents.data());
	}
	T.stop(); T.printTime("gid");
	
	cout << "Summary: \n";
	cout << "Particle pairs compared per mergeCell: " <<  float(particles_compared_per_cellpair)/count_pcpc << endl;
	cout << "Cells compared per cell: " <<  float(cells_compared_per_cell)/count_ccpc << endl;


}


void PointCloud::group_grid_hash2(float Rg){

	SimpleTimer T; T.reset(); T.start();
	calcGridParams(Rg/sqrt(3), par);
	T.stop(); T.printTime("calcGrid");

	
	float3 * pos = (float3*)points.data();

	T.reset(); T.start();	
	vector <unsigned int> point_hashes(nverts);
	vector <unsigned int> point_ids(nverts);
	// get the cell ID for each particle
	for (int i=0; i<nverts; ++i) {
		int3 cell_id = cellIndex(pos[i], par.origin);
		point_hashes[i] = cellHash(cell_id);
		point_ids[i]    = i;
	}
	T.stop(); T.printTime("cell Ids");


	// sort particles by cell ID
	T.reset(); T.start();	
	sort(point_ids.begin(), point_ids.end(), [&point_hashes](unsigned int i, unsigned int j){return point_hashes[i] < point_hashes[j];}); 
	T.stop(); T.printTime("sort");


	vector <int> parents(nverts);
	vector <int> sz(nverts,1);
	for (int i=0; i<nverts; ++i) parents[i] = i;


	map <int, int2> filled_cells;
	
	// Merge all points in a given cell into same group
	T.reset(); T.start();
	// start and end of each cell (both reflect indices in sorted hashes list*)
	//  -- *sorted hashes list is accessed as point_hashes[point_ids[0:nverts]]
	int start = 0, next=1;
	int cellStart = point_ids[start];
	while (next < nverts){
		int cellNext = point_ids[next];
		if (point_hashes[cellNext] != point_hashes[cellStart]){
			int2 se; se.x = start; se.y = next-1;
			filled_cells[point_hashes[cellStart]] = se;

			start = next;
			cellStart = cellNext; 
		}
		else{
			unite(cellNext, cellStart, parents.data(), sz.data());
		}
		++next;
	}
	T.stop(); T.printTime("cell s/e");


//	cout << "ALL:" << endl;
//	for (int i=0; i<30; ++i)
//	cout << point_ids[i] << " " << point_hashes[point_ids[i]] << " " << endl;

//	cout << "MAP (" << filled_cells.size() << "):" << endl;	
//	int i=0;
//	for (map<int,int2>::iterator it = filled_cells.begin(); it != filled_cells.end(); ++it){
//		cout << it->first << " " << it->second.x << " " << it->second.y << "\n";
//		++i; if (i >10) break;
//	}
//	cout << (--filled_cells.end())->first << " " << (--filled_cells.end())->second.x << " " << (--filled_cells.end())->second.y << endl;
//	// 

	T.reset(); T.start();
	
	int hashTable_size = 6000011;
	vector <HashNode> hashTable(hashTable_size);
	
//	vector <int> final_hash_ids(hashTable_size, -1);
	
//	vector <int3> cells;
//	cells.reserve(filled_cells.size());
	
	int avg_attempts = 0;
	int max_attempts = 0;
	for (map <int,int2>::iterator it = filled_cells.begin(); it != filled_cells.end(); ++it){
		
		int c1 = it->first;
		int3 cell;
		cell.z = int(c1/par.gridSize.x/par.gridSize.y);
		cell.y = int((c1 - cell.z*par.gridSize.x*par.gridSize.y)/par.gridSize.x);
		cell.x = int((c1 - cell.z*par.gridSize.x*par.gridSize.y - cell.y*par.gridSize.x));

		int2 val; val.x = it->second.x; val.y = it->second.y;
		int a = hash_insert(cell, val, hashTable.data(), hashTable_size);
		avg_attempts += a;
		max_attempts = max(max_attempts, a);
		
//		cells.push_back(cell);
	}
	cout << "Insertion Attempts: " << float(avg_attempts)/filled_cells.size() << " " << max_attempts << endl;
	T.stop(); T.printTime("map cells");
	

	
//	cout << "Counts:\n";
//	for (int i=0; i<hashTable_counts.size(); ++i){
//		if (hashTable_counts[i] == 0) cout << "- ";
//		else cout << hashTable_counts[i] << " ";	
//	}
//	cout << endl;
//	int b;
//	int avg_attempts2 = 0;
//	int max_attempts2 = 0;
////	cout << "Final IDs:\n";
//	for (int i=0; i<cells.size(); i=i+1){
//		int id_try = rand() % cells.size();
//		int id = hash_find(cells[id_try], hashTable.data(), hashTable_size, b);
////		cout << printCell(cells[id_try]) << ": " << final_hash_ids[id_try] << ": " 
////		     << printCell(hashTable[final_hash_ids[id_try]].key) << "  " <<  id << printCell(hashTable[id].key) << "\n";
//			avg_attempts2 += b;
//			max_attempts2 = max(max_attempts2, b);
//			if (id != final_hash_ids[id_try]) cout << "FATAL: Wrong ID Retreived" << endl;
//	}	
//	cout << "Retrieval Attempts (exisitng): " << float(avg_attempts2)/cells.size() << " " << max_attempts2 << endl;


//	int avg_attempts1 = 0;
//	int max_attempts1 = 0;
//	int attempts;
////	cout << "Random IDs:\n";
//	for (int i=0; i<1000000; ++i){
//		int3 cell;
//		cell.x = rand() % 1000;
//		cell.y = rand() % 1000;
//		cell.z = rand() % 1000;
//		int id = hash_find(cell, hashTable.data(), hashTable_size, attempts);
////		if (id >= 0) cout << printCell(cell) << ": " << "  " <<  id << printCell(hashTable_keys[id]) << "\n";
//		//else cout << "-" << " " << id << "\n";
//		avg_attempts1 += attempts;
//		max_attempts1 = max(max_attempts1, attempts);
//	}	
//	cout << "Retrieval Attempts (random): " << float(avg_attempts1)/1000000 << " " << max_attempts1 << endl;
	
	
	// group grid cells
	long long int pairs = 0;
	T.reset(); T.start();
	for (map<int,int2>::iterator it = filled_cells.begin(); it != filled_cells.end(); ++it){

		int c1 = it->first;

		int3 cell;
		cell.z = int(c1/par.gridSize.x/par.gridSize.y);
		cell.y = int((c1 - cell.z*par.gridSize.x*par.gridSize.y)/par.gridSize.x);
		cell.x = int((c1 - cell.z*par.gridSize.x*par.gridSize.y - cell.y*par.gridSize.x));

		int count=0;
		
		for (int ix = -2; ix <=2; ++ix){
			for (int iy = -2; iy <=2; ++iy){
				for (int iz = -2; iz <=2; ++iz){
					
					int3 cell_new;
					cell_new.x = clamp(cell.x + ix, 0, par.gridSize.x-1);
					cell_new.y = clamp(cell.y + iy, 0, par.gridSize.y-1);
					cell_new.z = clamp(cell.z + iz, 0, par.gridSize.z-1);
//					cout << "\tcell: " << cell_new.x << ", " << cell_new.y << ", " << cell_new.z << endl;

					int c2 = cellHash(cell_new);
					
					if (c2 < c1){	
//						int3 res = mergeCells_hash(cell, cell_new, Rg, hashTable.data(), point_ids, hashTable_size);
						int3 res = mergeCells(c1, c2, Rg, filled_cells, point_ids);
						if (res.x == 1) unite(res.y, res.z, parents.data(), sz.data());

						++count;
						++pairs;
					}
//					if (pairs % 10000 == 0) cout << pairs << "pairs compared.\n"; 
				}
			}
		}
		
		cells_compared_per_cell += count;
		count_ccpc += 1;
	}
	cout << "pairs = " << pairs << endl;
	T.stop(); T.printTime("group cells");
	
	T.reset(); T.start();
	group_ids.resize(nverts);
	for (int i=0; i<nverts; ++i){
		group_ids[i] = root(i, parents.data());
	}
	T.stop(); T.printTime("gid");
	
	cout << "Summary: \n";
	cout << "Particle pairs compared per mergeCell: " <<  float(particles_compared_per_cellpair)/count_pcpc << endl;
	cout << "Cells compared per cell: " <<  float(cells_compared_per_cell)/count_ccpc << endl;
	
	

}
	


void PointCloud::countNeighbours(int c1, int c2, float Rd, map <int, int2> & cells, vector <unsigned int> & pt_ids, vector <int> &n_nb){
	
	int count = 0;
	
	map <int, int2>::iterator ic1 = cells.find(c1);
	map <int, int2>::iterator ic2 = cells.find(c2);
	if (ic1 == cells.end() || ic2 == cells.end()) return; // if either cell doesnt exist, cant merge cells
	
	for (int p1 = ic1->second.x; p1 <= ic1->second.y; ++p1){	
		for (int p2 = ic2->second.x; p2 <= ic2->second.y; ++p2){
			++count;
			
			int ip1 = pt_ids[p1];
			int ip2 = pt_ids[p2];
			if (distance(ip1, ip2) < Rd){
				
				++n_nb[ip1];
				++n_nb[ip2];			
								
			} 
		}
	}
	
	particles_compared_per_cellpair += count;
	count_pcpc += 1;
	return;
}

int3 PointCloud::countNeighbours_hash(int3 c1, int3 c2, float Rd, HashNode * ht, vector <unsigned int> & pt_ids, int length, vector <int> &n_nb){
	
	int count = 0;
	
	int3 result; // (yes/no, point_id1, point_id2)
	result.x = result.y = result.z = 0;
	
	int attempts1, attempts2;
	int id_c1 = hash_find(c1, ht, length, &attempts1);
	int id_c2 = hash_find(c2, ht, length, &attempts2);
	if (id_c1 == -1 || id_c2 == -1) return result; // if either cell doesnt exist, cant merge cells
	
	for (int p1 = ht[id_c1].value.x; p1 <= ht[id_c1].value.y; ++p1){	
		for (int p2 = ht[id_c2].value.x; p2 <= ht[id_c2].value.y; ++p2){
			++count;
			
			int ip1 = pt_ids[p1];
			int ip2 = pt_ids[p2];
			if (distance(ip1, ip2) < Rd){
				++n_nb[ip1];
				++n_nb[ip2];			
			} 
		}
	}
//	cout << "attempts1 = " << attempts1 << " " << "attempts2 = " << attempts2 << endl;
	particles_compared_per_cellpair += count;
	count_pcpc += 1;
	return result;
}

void PointCloud::denoise(float Rd){

	SimpleTimer T; T.reset(); T.start();
	calcGridParams(Rd/sqrt(3), par);
	T.stop(); T.printTime("calcGrid");

	
	float3 * pos = (float3*)points.data();

	T.reset(); T.start();	
	vector <unsigned int> point_hashes(nverts);
	vector <unsigned int> point_ids(nverts);
	// get the cell ID for each particle
	for (int i=0; i<nverts; ++i) {
		int3 cell_id = cellIndex(pos[i], par.origin);
		point_hashes[i] = cellHash(cell_id);
		point_ids[i]    = i;
	}
	T.stop(); T.printTime("cell Ids");


	// sort particles by cell ID
	T.reset(); T.start();	
	sort(point_ids.begin(), point_ids.end(), [&point_hashes](unsigned int i, unsigned int j){return point_hashes[i] < point_hashes[j];}); 
	T.stop(); T.printTime("sort");

	vector <int> parents(nverts);
	vector <int> sz(nverts,1);
	for (int i=0; i<nverts; ++i) parents[i] = i;

//	cout << "ALL:" << endl;
//	for (int i=0; i<30; ++i)
//	cout << point_ids[i] << " " << point_hashes[point_ids[i]] << " " << endl;


	int hashTable_size = 6000011;
	vector <HashNode> hashTable(hashTable_size);

	int n_attempts = 0;	
	int avg_attempts = 0;
	int max_attempts = 0;

	// Merge all points in a given cell into same group
	T.reset(); T.start();
	// start and end of each cell (both reflect indices in sorted hashes list*)
	//  -- *sorted hashes list is accessed as point_hashes[point_ids[0:nverts]]
	int start = 0, next=1;
	int cellStart = point_ids[start];
	while (next < nverts){
		int cellNext = point_ids[next];
		if (point_hashes[cellNext] != point_hashes[cellStart]){
			int2 se; se.x = start; se.y = next-1;
			
			int3 cell;
			cell.z = int(point_hashes[cellStart]/par.gridSize.x/par.gridSize.y);
			cell.y = int((point_hashes[cellStart] - cell.z*par.gridSize.x*par.gridSize.y)/par.gridSize.x);
			cell.x = int((point_hashes[cellStart] - cell.z*par.gridSize.x*par.gridSize.y - cell.y*par.gridSize.x));

			int a = hash_insert(cell, se, hashTable.data(), hashTable_size);
			++n_attempts;
			avg_attempts += a;
			max_attempts = max(max_attempts, a);

			start = next;
			cellStart = cellNext; 
		}
		else {
			//unite(cellNext, cellStart, parents.data(), sz.data());
			//cout << "uniting " << cellStart << "|" << cellNext << " (" << point_hashes[cellStart] << "|" << point_hashes[cellNext] << endl;  
		}
		++next;
	}
	cout << "Insertion Attempts: " << float(avg_attempts)/n_attempts << " " << max_attempts << endl;
	T.stop(); T.printTime("map cells");


//	cout << "ALL:" << endl;
//	for (int i=0; i<30; ++i)
//	cout << point_ids[i] << " " << point_hashes[point_ids[i]] << " " << group_ids[point_ids[i]] << endl;

//	cout << "MAP (" << filled_cells.size() << "):" << endl;	
//	int i=0;
//	for (map<int,int2>::iterator it = filled_cells.begin(); it != filled_cells.end(); ++it){
//		cout << it->first << " " << it->second.x << " " << it->second.y << "\n";
//		++i; if (i >10) break;
//	}
//	cout << (--filled_cells.end())->first << " " << (--filled_cells.end())->second.x << " " << (--filled_cells.end())->second.y << endl;
//	// 

	neighbourCounts.resize(nverts,0);

	particles_compared_per_cellpair = 0;
	count_pcpc = 0;

	// group grid cells
	int n_filled_cells = 0;
	long long int pairs = 0;
	T.reset(); T.start();
	for (int i=0; i<hashTable_size; ++i){
		
		int3 cell = hashTable[i].key;
		if (cell == HashNode().key) continue;
		++n_filled_cells;
		
		int count=0;
		
		for (int ix = -2; ix <=2; ++ix){
			for (int iy = -2; iy <=2; ++iy){
				for (int iz = -2; iz <=2; ++iz){
					
					int3 cell_new;
					cell_new.x = clamp(cell.x + ix, 0, par.gridSize.x-1);
					cell_new.y = clamp(cell.y + iy, 0, par.gridSize.y-1);
					cell_new.z = clamp(cell.z + iz, 0, par.gridSize.z-1);
//					cout << "\tcell: " << cell_new.x << ", " << cell_new.y << ", " << cell_new.z << endl;

					if (cellHash(cell_new) <= cellHash(cell)){	
//						int3 res = mergeCells(c1, c2, Rg, filled_cells, point_ids);
						countNeighbours_hash(cell, cell_new, Rd, hashTable.data(), point_ids, hashTable_size, neighbourCounts);

						++count;
						++pairs;
					}
//					if (pairs % 10000 == 0) cout << pairs << "pairs compared.\n"; 
				}
			}
		}
		
		cells_compared_per_cell += count;
		count_ccpc += 1;
	}
	cout << "pairs = " << pairs << endl;
	T.stop(); T.printTime("group cells");
	
	T.reset(); T.start();
	group_ids.resize(nverts);
	for (int i=0; i<nverts; ++i){
		group_ids[i] = root(i, parents.data());
	}
	T.stop(); T.printTime("gid");
	
	cout << "Summary: \n";
	cout << "Particle pairs compared per mergeCell: " <<  float(particles_compared_per_cellpair)/count_pcpc << endl;
	cout << "Cells compared per cell: " <<  float(cells_compared_per_cell)/count_ccpc << endl;
	cout << "No. of filled cells = " << n_filled_cells << endl;
	

}





