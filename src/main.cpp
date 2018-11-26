#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <cmath>
#include <liblas/liblas.hpp>
using namespace std;

#include "../headers/graphics.h"
#include "../utils/simple_math.h"

namespace fl{

struct vec3{
	float x;
	float y;
	float z;
};

}


struct float3{
	float x;
	float y;
	float z;
};

struct int3{
	int x;
	int y;
	int z;
};

struct int2{
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
	
	void deleteGround(){
		for (int k=0; k<nverts; ++k){
			int ix = floor((points[3*k+0]- dem.xmin)/ dem.dx);
			int iy = floor((points[3*k+1]- dem.ymin)/ dem.dy);
			if (points[3*k+2] < 0.5) points[3*k+0] = points[3*k+1] = points[3*k+2] = 0; 
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

};


int main(int argc, char **argv){

	PointCloud cr;
//	cr.read_las("/home/jaideep/Data/LIDAR/2017-02-20_21-47-24.las");
//	cr.createDEM(0.5,0.5);
//	cr.dem.printToFile("dem.txt");
//	cr.subtractDEM();
//	cr.deleteGround();
	cr.generateRandomClusters(4000000, 500, -100, 100, -100, 100, 0, 10, 2);
	
	init_hyperGL(&argc, argv);

	Palette p(100);
	p.create_rainbow();
	p.create_random();
	
	cout << "sort...\n";
	sort_by_z(cr.points.data(), cr.points.data()+3*cr.nverts); 
	for (int i=0; i<10; ++i) cout << cr.points[3*i] << " " << cr.points[3*i+1] << " " << cr.points[3*i+2] << endl;

//	cr.group_serial(2);
	cr.group_grid(0.5);

//	vector <float> cols9z = p.map_values(&cr.points[2], cr.nverts, 3);	// map z value
	vector <float> gids(cr.group_ids.begin(), cr.group_ids.end());
	vector <float> cols9z = p.map_values(gids.data(), cr.nverts, 1);	// map group ID
	Shape pt(cr.nverts, 3, "points", true); //, 4, -1, 1);
	pt.setVertices(cr.points.data());	
	pt.setColors(&cols9z[0]);
	vector <float> ex = calcExtent(cr.points.data(), cr.nverts, 3);
	pt.setExtent(ex);
//	pt.pointSize = 2;
 
 
//	vector <int> slices = z_slices(cr.points.data(), cr.nverts, 0.1);
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




void PointCloud::group_hashedGrid(float Rg){
	
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
		point_hashes[i] = calcHash(cell_id);
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
	vector <int> cell_hashes;
	vector <int> cell_start;
	vector <int> cell_end;
	cell_hashes.reserve(nverts);	// number of non-empty cells can never exceed nverts
	cell_start.reserve(nverts);		// number of non-empty cells can never exceed nverts
	cell_end.reserve(nverts);		// number of non-empty cells can never exceed nverts
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
					
					unsigned int hash = calcHash(cell_new);
					
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





