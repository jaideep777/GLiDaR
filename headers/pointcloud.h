#ifndef POINT_CLOUD_H
#define POINT_CLOUD_H
#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include "hashtable.h"

using namespace std;

#ifdef __CUDACC__

#define DEVICE_NAME __device__ __host__

#else

#define DEVICE_NAME 

class float3{
	public:
	float x;
	float y;
	float z;
};


#endif

namespace fl{

class vec3{
	public:
	float x;
	float y;
	float z;
};

}



struct Params{
	float3 cellSize;
	int3 gridSize;
	float3 origin;
};


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
	
	void deleteGround(float ht);
	
	void generateRandomClusters(int N, int nc, float xmin, float xmax, float ymin, float ymax, float zmin, float zmax, float R);


	public:
	vector<int> group_ids; 
	float distance(int p, int q);
	
	void group_serial(float Rg);
	
	void calcGridParams(float Rg, Params & par);
	void group_grid(float Rg);
	
	int3 mergeCells(int c1, int c2, float Rg, map <int, int2> & cells, vector <unsigned int> & pt_ids);	
	void group_grid_map(float Rg);

	int3 mergeCells_hash(int3 c1, int3 c2, float Rg, HashNode * ht, vector <unsigned int> & pt_ids, int length);
	void group_grid_hash(float Rg);
	void group_grid_hash2(float Rg);

	int3 mergeCells_hash_stl(int3 c1, int3 c2, float Rg, unordered_map<int3,int2,KeyHasher> &ht, vector <unsigned int> & pt_ids, int length);
	void group_grid_hashSTL(float Rg);


	vector <int> noisePoints;
	void countNeighbours(int c1, int c2, float Rd, map <int, int2> & cells, vector <unsigned int> & pt_ids, vector <int> &n_nb);
	void denoise(float Rd);

};



#endif


