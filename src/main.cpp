#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <cmath>
#include <liblas/liblas.hpp>
using namespace std;

#include "../headers/graphics.h"
#include "../utils/simple_math.h"

namespace fl{

class vec3{
	public:
	float x;
	float y;
	float z;
};

}


class float3{
	public:
	float x;
	float y;
	float z;
};

class int3{
	public:
	int x;
	int y;
	int z;
	
	bool operator != (int3 v){
		return (x != v.x || y != v.y || z != v.z);
	}
};

class int2{
	public:
	int x, y;
};


struct Params{
	float3 cellSize;
	int3 gridSize;
	float3 origin;
} par;


bool compare_z(const fl::vec3& p1, const fl::vec3& p2){
	return (p1.z < p2.z);
}

bool compare_y(const fl::vec3& p1, const fl::vec3& p2){
	return (p1.y < p2.y);
}

bool compare_x(const fl::vec3& p1, const fl::vec3& p2){
	return (p1.x < p2.x);
}

void sort_by_y(float* begin, float* end){
	sort((fl::vec3*)begin, (fl::vec3*)end, compare_y); 
}
	
void sort_by_z(float* begin, float* end){
	sort((fl::vec3*)begin, (fl::vec3*)end, compare_z); 
}

void sort_by_x(float* begin, float* end){
	sort((fl::vec3*)begin, (fl::vec3*)end, compare_x); 
}


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











class DigitalElevModel{
	public:
	int   nx, ny;	// number of rows and columns in DEM matrix
	float dx, dy;	// x and y resolution (size of each grid-cell in DEM matrix)
	float x0, y0;	// centre of the lower-left gridcell
	float xmin, xmax, ymin, ymax;	// bounds
	vector<float>elev;				// elevation

	public:
//	DigitalElevModel(int _nx, int _ny, float _dx, float _dy, float _x0, float _y0){
//		nx   = _nx;			
//		ny   = _ny;			
//		dx   = _dx;			
//		dy   = _dy;			
//		x0   = _x0;			
//		y0   = _y0;	
//			
//	}

	DigitalElevModel(){}

	void init(int _nx, int _ny, float _dx, float _dy, float _x0, float _y0){
		nx   = _nx;			
		ny   = _ny;			
		dx   = _dx;			
		dy   = _dy;			
		x0   = _x0;			
		y0   = _y0;
		
		xmin = x0 - dx/2;
		ymin = x0 - dx/2;
		xmax = (dx*nx) + xmin; 
		ymax = (dy*ny) + ymin;
		
		elev.resize(nx*ny);
		
	}
	
	void initFromBounds(int _xmin, int _xmax, float _ymin, float _ymax, float _dx, float _dy){
		xmin = _xmin;			
		xmax = _xmax;			
		ymin = _ymin;			
		ymax = _ymax;			
		dx   = _dx;			
		dy   = _dy;
		nx   =  (xmax- xmin)/ dx ; 
		ny   =  (ymax- ymin)/ dy ; 
		x0   =  xmin + dx/2; 
		y0   =  ymin + dy/2;
		
		elev.resize(nx*ny);
		 			
	}
	
	vector<int> get_index(float x, float y){
		vector<int>index(2);
		index[0] = floor((x- xmin)/ dx);
		index[1] = floor((y- ymin)/ dy);
		return index;
	}
	
	void printToFile(string filename){
		ofstream fout(filename.c_str());
		for (int j=0; j<ny; ++j){
			for (int i=0; i<nx; ++i){
				fout << elev[j*nx+i] << "\t";
			}
			fout << "\n";
		}
		fout.close();
	}
	
};


class PointCloud{
	public:
	int nverts;
//	float dx,dy;
//	int nx,ny;
	vector <float> points;

	DigitalElevModel dem;

	public:
	void read_las(string file){
		ifstream ifs;
		ifs.open(file.c_str(),ios::in | ios::binary);
		if (!ifs) cout << "Error Opening file: " << file << endl;

		liblas::ReaderFactory f;
		liblas::Reader reader = f.CreateWithStream(ifs);
		liblas::Header const& header = reader.GetHeader();
		cout << "Compressed: " << ((header.Compressed() == true) ? "true\n":"false\n");
		cout << "Signature: " << header.GetFileSignature() << '\n';
		cout << "Points count: " << header.GetPointRecordsCount() << '\n';
		nverts = header.GetPointRecordsCount();

		points.reserve(3*nverts);
		
		while (reader.ReadNextPoint()){ //count<10
			liblas::Point const& p = reader.GetPoint();
			points.push_back(p.GetX());
			points.push_back(p.GetY());
			points.push_back(p.GetZ());

			//cout << p.GetX() << ", " << p.GetY() << ", " << p.GetZ() << "\n";
		}

	}
	
	void createDEM(float dx, float dy){
		float xmin = min_element((fl::vec3*)points.data(), (fl::vec3*)(points.data()+3*nverts), compare_x)->x;
		float xmax = max_element((fl::vec3*)points.data(), (fl::vec3*)(points.data()+3*nverts), compare_x)->x;
		float ymin = min_element((fl::vec3*)points.data(), (fl::vec3*)(points.data()+3*nverts), compare_y)->y;
		float ymax = max_element((fl::vec3*)points.data(), (fl::vec3*)(points.data()+3*nverts), compare_y)->y;
		cout << "x range: " << xmin << " " << xmax << endl;	  
		cout << "y range: " << ymin << " " << ymax << endl;		
		
		xmin = floor(xmin/dx)*dx;
		xmax = ceil(xmax/dx)*dx;
		ymin = floor(ymin/dy)*dy;
		ymax = ceil(ymax/dy)*dy;

		dem.initFromBounds(xmin,xmax,ymin,ymax,dx,dy);
		//dem.init(nx,ny,dx,dy,xmin,ymin);
		
		fill(dem.elev.begin(), dem.elev.end(), 1e20);	// initialize all elevations with 1e20
		for (int k=0; k<nverts; ++k){
			int ix = floor((points[3*k+0]- dem.xmin)/ dem.dx);
			int iy = floor((points[3*k+1]- dem.ymin)/ dem.dy);
			dem.elev[iy*dem.nx+ix] = min(points[3*k+2], dem.elev[iy*dem.nx+ix]);
		}
		
		// TODO: Filter DEM
	}

	void subtractDEM(){
		for (int k=0; k<nverts; ++k){
			int ix = floor((points[3*k+0]- dem.xmin)/ dem.dx);
			int iy = floor((points[3*k+1]- dem.ymin)/ dem.dy);
			points[3*k+2] -= dem.elev[iy*dem.nx+ix];	// TODO: subract interpolated value
		}
	}	
	
	void deleteGround(float ht){
		for (int k=0; k<nverts; ++k){
			int ix = floor((points[3*k+0]- dem.xmin)/ dem.dx);
			int iy = floor((points[3*k+1]- dem.ymin)/ dem.dy);
			if (points[3*k+2] < ht) points[3*k+0] = points[3*k+1] = points[3*k+2] = 0; 
		}
	}
	
	void generateRandomClusters(int N, int nc, float xmin, float xmax, float ymin, float ymax, float zmin, float zmax, float R){
		nverts = N;
		points.resize(3*nverts);
		
		// create centroids
		vector <float> centroids;
		centroids.reserve(nc);
		for (int i=0; i<nc; ++i){
			centroids.push_back(runif(xmin, xmax));
			centroids.push_back(runif(ymin, ymax));
			centroids.push_back(runif(zmin, zmax));
		}

		for (int i=0; i<N; ++i){
			int g = rand() % nc;
			float dx = runif(-R, R);
			float dy = runif(-R, R);
			float dz = runif(-R, R);
		
			points[3*i+0] = centroids[3*g+0] + dx;
			points[3*i+1] = centroids[3*g+1] + dy;
			points[3*i+2] = centroids[3*g+2] + dz;
		}
	}



	vector<int> group_ids; 
	float distance(int p, int q);
	
	void group_serial(float Rg);
	
	void calcGridParams(float Rg, Params & par);
	void group_grid(float Rg);
	
	int3 mergeCells(int c1, int c2, float Rg, map <int, int2> & cells, vector <unsigned int> & pt_ids);	
	void group_grid_map(float Rg);
	void group_grid_hash(float Rg);

	vector <int> noisePoints;
	void countNeighbours(int c1, int c2, float Rd, map <int, int2> & cells, vector <unsigned int> & pt_ids, vector <int> &n_nb);
	void denoise(float Rd);

};


int main(int argc, char **argv){

	PointCloud cr;
	cr.read_las("/home/jaideep/Data/LIDAR/2017-02-20_21-47-24.las");
	cr.createDEM(0.5,0.5);
	cr.dem.printToFile("dem.txt");
	cr.subtractDEM();
//	cr.deleteGround(0.2);
//	cr.generateRandomClusters(1000000, 500, -100, 100, -100, 100, 0, 10, 2);
	
	init_hyperGL(&argc, argv);

	Palette p(100);
	p.create_rainbow();
//	p.create_random();
	
	cout << "sort...\n";
	sort_by_z(cr.points.data(), cr.points.data()+3*cr.nverts); 
	for (int i=0; i<10; ++i) cout << cr.points[3*i] << " " << cr.points[3*i+1] << " " << cr.points[3*i+2] << endl;

//	cr.group_serial(2);
//	cr.group_grid_map(2);
//	cr.denoise(2);
	cr.group_grid_hash(0.1);

	vector <float> cols9z = p.map_values(&cr.points[2], cr.nverts, 3);	// map z value
//	vector <float> gids(cr.group_ids.begin(), cr.group_ids.end());
//	vector <float> cols9z = p.map_values(gids.data(), cr.nverts, 1);	// map group ID
	Shape pt(cr.nverts, 3, "points", true); //, 4, -1, 1);
	pt.setVertices(cr.points.data());	
	pt.setColors(&cols9z[0]);
	vector <float> ex = calcExtent(cr.points.data(), cr.nverts, 3);
	pt.setExtent(ex);
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
	
	Shape axis(6, 3, "lines");
	axis.setVertices(pos3);
	axis.setColors((float*)col3);
	
	
	while(1){	   // infinite loop needed to poll anim_on signal.
//		slices_shapes[current_render]->b_render = false;
//		current_render = generic_count % (slices_shapes.size()-1);
////		cout << "Currently rendering: " << current_render << "\n";
//		slices_shapes[current_render]->b_render = true;
		
		glutMainLoopEvent();
		usleep(20000);
	}
	// launch sim end.
	
	return 0;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Union-Find functions to calculate group Indices 
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// This file contains 3 versions of union-find:
// 		1. serial version is the standard UF algorithm. Note that latest positions 
//		   must be copied to host before this function can be called.
// 		2. Parallel version with pairwise comparisons in parallel followed by  
//		   serial looping over all pairs
//		3. Parallel_sort version with parallel pairwise comparisons followed
//		   by atomic sorting of pairs, then serial unites over only the 
//		   close pairs. This is the fastest version.
//
// Serial version supports constant and variable grouping radius 
// Parallel versions take a constant radius of grouping from Movementparams. Variable
// 		grouping radius code is not implemented in the parallel version because
// 		it is not used anyway for the current simulations.
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

float PointCloud::distance(int p, int q){
	float dx = points[3*p+0] - points[3*q+0];
	float dy = points[3*p+1] - points[3*q+1];
	float dz = points[3*p+2] - points[3*q+2];
	return sqrt(dx*dx + dy*dy + dz*dz);
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







void PointCloud::calcGridParams(float Rg, Params & par){
	float3 * pos = (float3*)points.data();
	par.cellSize.x = par.cellSize.y = par.cellSize.z = Rg;
	
	float xmin = min_element((fl::vec3*)pos, (fl::vec3*)(pos+nverts), compare_x)->x;
	float xmax = max_element((fl::vec3*)pos, (fl::vec3*)(pos+nverts), compare_x)->x;
	float ymin = min_element((fl::vec3*)pos, (fl::vec3*)(pos+nverts), compare_y)->y;
	float ymax = max_element((fl::vec3*)pos, (fl::vec3*)(pos+nverts), compare_y)->y;
	float zmin = min_element((fl::vec3*)pos, (fl::vec3*)(pos+nverts), compare_z)->z;
	float zmax = max_element((fl::vec3*)pos, (fl::vec3*)(pos+nverts), compare_z)->z;

	cout << "Group: " << endl;
	cout << "x range: " << xmin << " " << xmax << endl;	  
	cout << "y range: " << ymin << " " << ymax << endl;		
	cout << "z range: " << zmin << " " << zmax << endl;		
	
	xmin = floor(xmin/Rg)*Rg;
	xmax = ceil(xmax/Rg)*Rg;
	ymin = floor(ymin/Rg)*Rg;
	ymax = ceil(ymax/Rg)*Rg;
	zmin = floor(zmin/Rg)*Rg;
	zmax = ceil(zmax/Rg)*Rg;

	float3 origin;
	origin.x = xmin; origin.y = ymin; origin.z = zmin;
	par.origin = origin;

	par.gridSize.x  =  (xmax- xmin)/ Rg ; 
	par.gridSize.y  =  (ymax- ymin)/ Rg ; 
	par.gridSize.z  =  (zmax- zmin)/ Rg ; 

	cout << "Group: " << endl;
	cout << "x range: " << xmin << " " << xmax << endl;	  
	cout << "y range: " << ymin << " " << ymax << endl;		
	cout << "z range: " << zmin << " " << zmax << endl;		

	cout << "gridSize = " << par.gridSize.x*par.gridSize.y*par.gridSize.z << endl;
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
					cell_new.z = clamp(cell.z + iz, 0, par.gridSize.z-1);
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


void PointCloud::denoise(float Rd){

	SimpleTimer T; T.reset(); T.start();
	calcGridParams(Rd, par);
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
		++next;
	}
	T.stop(); T.printTime("cell s/e");


	cout << "ALL:" << endl;
	for (int i=0; i<30; ++i)
	cout << point_ids[i] << " " << point_hashes[point_ids[i]] << " " << endl;

	cout << "MAP (" << filled_cells.size() << "):" << endl;	
	int i=0;
	for (map<int,int2>::iterator it = filled_cells.begin(); it != filled_cells.end(); ++it){
		cout << it->first << " " << it->second.x << " " << it->second.y << "\n";
		++i; if (i >10) break;
	}
	cout << (--filled_cells.end())->first << " " << (--filled_cells.end())->second.x << " " << (--filled_cells.end())->second.y << endl;
	// 

	noisePoints.resize(nverts,0);
	
//	// find cells with no neighbours
	long long int pairs = 0;
	T.reset(); T.start();
	for (map<int,int2>::iterator it = filled_cells.begin(); it != filled_cells.end(); ++it){

		int c1 = it->first;

		int3 cell;
		cell.z = int(c1/par.gridSize.x/par.gridSize.y);
		cell.y = int((c1 - cell.z*par.gridSize.x*par.gridSize.y)/par.gridSize.x);
		cell.x = int((c1 - cell.z*par.gridSize.x*par.gridSize.y - cell.y*par.gridSize.x));

		int count=0;
		
		for (int ix = -1; ix <=1; ++ix){
			for (int iy = -1; iy <=1; ++iy){
				for (int iz = -1; iz <=1; ++iz){
					
					int3 cell_new;
					cell_new.x = clamp(cell.x + ix, 0, par.gridSize.x-1);
					cell_new.y = clamp(cell.y + iy, 0, par.gridSize.y-1);
					cell_new.z = clamp(cell.z + iz, 0, par.gridSize.z-1);
//					cout << "\tcell: " << cell_new.x << ", " << cell_new.y << ", " << cell_new.z << endl;

					int c2 = cellHash(cell_new);
					
					if (c2 < c1){	
						countNeighbours(c1, c2, Rd, filled_cells, point_ids, noisePoints);

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
	
	
	cout << "Summary: \n";
	cout << "Particle pairs compared per mergeCell: " <<  float(particles_compared_per_cellpair)/count_pcpc << endl;
	cout << "Cells compared per cell: " <<  float(cells_compared_per_cell)/count_ccpc << endl;


}


// HASH TABLE IMPLEMENTATION OF GRID

int hash3D(int3 key, int length){
	return (key.z*377771 + key.y*677 + key.x) % length;	// 2.15, 64
//	return (key.z*377771 + key.y*133337 + key.x) % length; // 1.7, 25
//	return (key.z*497771 + key.y*133337 + key.x) % length; // 30, 173
//	return ((long long int)key.z*497771 + key.y*787 + key.x) % length; // 1.5, 28
}

int insert(int3 key, int value, int3 * keys, int* values, int* counts, int length, int* final_id = NULL){
	int id = hash3D(key, length);
	int hash = id;
	int count = 1;
	while (values[id] != -1){ // while counts[id] != 0
		id = (id+11)%length;
		++count;
	}
	values[id] = value;
	keys[id] = key;
	counts[hash] = max(counts[id], count);
	*final_id = id;
	return count;

}


int hash_find(int3 key, int3 * keys, int* values, int* counts, int length, int& attempts){
	int id = hash3D(key, length);
	int hash = id;
	int count = 1; attempts = count;
	while (keys[id] != key){ //(keys[id].x != key.x || keys[id].y != key.y || keys[id].z != key.z){
		id = (id+11)%length;
		++count;
		attempts = count;
		if (count > counts[hash]) return -1;
	}
	return id;
}

#define printCell(cell) "(" << cell.x << "," << cell.y << "," << cell.z << ")"


void PointCloud::group_grid_hash(float Rg){

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


	// sort particles by cell ID
	T.reset(); T.start();	
	sort(point_ids.begin(), point_ids.end(), [&point_hashes](unsigned int i, unsigned int j){return point_hashes[i] < point_hashes[j];}); 
	T.stop(); T.printTime("sort");


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
		++next;
	}
	T.stop(); T.printTime("cell s/e");


	cout << "ALL:" << endl;
	for (int i=0; i<30; ++i)
	cout << point_ids[i] << " " << point_hashes[point_ids[i]] << " " << endl;

	cout << "MAP (" << filled_cells.size() << "):" << endl;	
	int i=0;
	for (map<int,int2>::iterator it = filled_cells.begin(); it != filled_cells.end(); ++it){
		cout << it->first << " " << it->second.x << " " << it->second.y << "\n";
		++i; if (i >10) break;
	}
	cout << (--filled_cells.end())->first << " " << (--filled_cells.end())->second.x << " " << (--filled_cells.end())->second.y << endl;
	// 

	T.reset(); T.start();
	vector <int3> cells;
	cells.reserve(filled_cells.size());
	
	for (map <int,int2>::iterator it = filled_cells.begin(); it != filled_cells.end(); ++it){
		
		int c1 = it->first;
		int3 cell;
		cell.z = int(c1/par.gridSize.x/par.gridSize.y);
		cell.y = int((c1 - cell.z*par.gridSize.x*par.gridSize.y)/par.gridSize.x);
		cell.x = int((c1 - cell.z*par.gridSize.x*par.gridSize.y - cell.y*par.gridSize.x));

		cells.push_back(cell);
	}
	T.stop(); T.printTime("map cells");
	
	
	cout << "CELLS VEC:\n"; 
	for (int i=0; i<10; ++i){
		cout << cellHash(cells[i]) << endl;
	}

	int hashTable_size = 6000011;
	vector <int3> hashTable_keys(hashTable_size);
	vector <int>  hashTable_vals(hashTable_size, -1);
	vector <int>  hashTable_counts(hashTable_size, 0);
	
	vector <int> final_hash_ids(hashTable_size, -1);
	
	int avg_attempts = 0;
	int max_attempts = 0;
	for (int i=0; i<cells.size(); ++i){
		int a = insert(cells[i], 1, hashTable_keys.data(), hashTable_vals.data(), hashTable_counts.data(), hashTable_size, &final_hash_ids[i]);
		avg_attempts += a;
		max_attempts = max(max_attempts, a);
	}
	cout << "Insertion Attempts: " << float(avg_attempts)/cells.size() << " " << max_attempts << endl;
	
//	cout << "Counts:\n";
//	for (int i=0; i<hashTable_counts.size(); ++i){
//		if (hashTable_counts[i] == 0) cout << "- ";
//		else cout << hashTable_counts[i] << " ";	
//	}
//	cout << endl;
	int b;
	int avg_attempts2 = 0;
	int max_attempts2 = 0;
//	cout << "Final IDs:\n";
	for (int i=0; i<cells.size(); i=i+1){
		int id = hash_find(cells[rand() % cells.size()], hashTable_keys.data(), hashTable_vals.data(), hashTable_counts.data(), hashTable_size, b);
//		cout << printCell(cells[i]) << ": " << final_hash_ids[i] << ": " 
//		     << printCell(hashTable_keys[final_hash_ids[i]]) << "  " <<  id << printCell(hashTable_keys[id]) << "\n";
			avg_attempts2 += b;
			max_attempts2 = max(max_attempts2, b);
	}	
	cout << "Retrieval Attempts (exisitng): " << float(avg_attempts2)/cells.size() << " " << max_attempts2 << endl;


	int avg_attempts1 = 0;
	int max_attempts1 = 0;
	int attempts;
//	cout << "Random IDs:\n";
	for (int i=0; i<1000000; ++i){
		int3 cell;
		cell.x = rand() % 1000;
		cell.y = rand() % 1000;
		cell.z = rand() % 1000;
		int id = hash_find(cell, hashTable_keys.data(), hashTable_vals.data(), hashTable_counts.data(), hashTable_size, attempts);
//		if (id >= 0) cout << printCell(cell) << ": " << "  " <<  id << printCell(hashTable_keys[id]) << "\n";
		//else cout << "-" << " " << id << "\n";
		avg_attempts1 += attempts;
		max_attempts1 = max(max_attempts1, attempts);
	}	
	cout << "Retrieval Attempts (random): " << float(avg_attempts1)/1000000 << " " << max_attempts1 << endl;
	

}
	
//int hash3D(int3 key, int length){
//	return (key.z*377771 + key.y*677 + key.x) % length;
//}

//int insert(int3 key, int value, int3 * keys, int* values, int length){
//	int id = hash3D(key, length);
//	int count = 1;
//	while (values[id] != -1){
//		id = (id+1)%length;
//		++count;
//	}
//	values[id] = value;
//	keys[id] = key;
//	return count;
//}

//int find(int3 key, int3 * keys, int* values, int length){
//	int id = hash3D(key, length);
//	while (keys[id].x != key.x && keys[id].y != key.y && keys[id].z != key.z){
//		id = (id+1)%length;
//	}
//	return values[id];
//}





