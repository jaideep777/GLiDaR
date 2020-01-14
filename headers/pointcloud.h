#ifndef POINT_CLOUD_H
#define POINT_CLOUD_H
#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include "hashtable.h"

using namespace std;

#include "../utils/simple_math.h"
#include "grid.h"



namespace fl{

class vec3{
	public:
	float x;
	float y;
	float z;
};

}

bool compare_z(const fl::vec3& p1, const fl::vec3& p2);
bool compare_y(const fl::vec3& p1, const fl::vec3& p2);
bool compare_x(const fl::vec3& p1, const fl::vec3& p2);

void sort_by_z(float* begin, float* end);
void sort_by_y(float* begin, float* end);	
void sort_by_x(float* begin, float* end);



class DigitalElevModel{
	public:
	int   nx, ny;	// number of rows and columns in DEM matrix
	float dx, dy;	// x and y resolution (size of each grid-cell in DEM matrix)
	float x0, y0;	// centre of the lower-left gridcell
	float xmin, xmax, ymin, ymax;	// bounds
	vector<float>elev;				// elevation

	public:

//	DigitalElevModel();

	void init(int _nx, int _ny, float _dx, float _dy, float _x0, float _y0);
	
	void initFromBounds(int _xmin, int _xmax, float _ymin, float _ymax, float _dx, float _dy);
	
	vector<int> get_index(float x, float y);
	
	void printToFile(string filename);
	
};



class PointCloud{
	public:
	int nverts;
//	float dx,dy;
//	int nx,ny;
	vector <float> points;

	DigitalElevModel dem;

	public:
	void read_las(string file);
	
	void createDEM(float dx, float dy);

	void subtractDEM();
	void subtractDEM_bil();
	
	void deleteGround(float ht);
	
	void generateRandomClusters(int N, int nc, float xmin, float xmax, float ymin, float ymax, float zmin, float zmax, float R);


	public:
	vector<int> group_ids; 
	float distance(int p, int q);
	
	void group_serial(float Rg);
	
	void calcGridParams(float Rg, Grid& par);
	void group_grid(float Rg);
	
	int3 mergeCells(int c1, int c2, float Rg, map <int, int2> & cells, vector <unsigned int> & pt_ids);	
	void group_grid_map(float Rg);

	int3 mergeCells_hash(int3 c1, int3 c2, float Rg, HashNode * ht, vector <unsigned int> & pt_ids, int length);
	void group_grid_hash(float Rg);
	void group_grid_hash2(float Rg);

	void group_grid_hash_gpu(float Rg);

	int3 mergeCells_hash_stl(int3 c1, int3 c2, float Rg, unordered_map<int3,int2,KeyHasher> &ht, vector <unsigned int> & pt_ids, int length);
	void group_grid_hashSTL(float Rg);


	vector <int> neighbourCounts;
	void countNeighbours(int c1, int c2, float Rd, map <int, int2> & cells, vector <unsigned int> & pt_ids, vector <int> &n_nb);
	int3 countNeighbours_hash(int3 c1, int3 c2, float Rd, HashNode * ht, vector <unsigned int> & pt_ids, int length, vector <int> &n_n);
	void denoise(float Rd);
	
	void countNeighbours_hash_gpu(float Rd);

};



#endif


