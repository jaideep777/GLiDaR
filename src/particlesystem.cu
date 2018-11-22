#include "particlesystem.h"
using namespace std;


struct par{
	float3 cellSize;
	int3 gridSize;	// if grid size is same as grouping distance, then only particles in same or neighnouring cells will qualify for union
}


// Get the 3D index of the cell that contains the particle
__device__ int3 getIndex(float3 pos, float3 origin){
	int3 cell;
	cell.x = (pos.x - origin.x)/par.cellSize.x;
	cell.y = (pos.y - origin.y)/par.cellSize.y;
	cell.z = (pos.z - origin.z)/par.cellSize.z;
	return cell;
}

// calculate the hash of a given cell
__device__ int3 cellHash(int3 cell){
	//		   iz * nx * ny                     +    iy * nx            + ix
	return cell.z*par.gridSize.x*par.gridSize.y + cell.y*par.gridSize.x + cell.x;
	// can use z-order curve as hash here rather than 1D cell-index
}

// calculate the hashes of all particles (i.e. hashes of cells containing those particles)
// store the [particle id, hash] pairs (in separate arrays, index of array gives the pair id)  
__global__ calcHash(float3 * pos, uint * hashes, uint * pids, int N){
	
	uint tid = blockIdx.x*blockDim.x + threadIdx.x;
	if (tid >= N) return;
	
	int3 cell = getIndex(pos[tid]);
	uint hash = cellHash(cell);
	
	hashes[tid] = hash;
	pids[tid]   = tid;
	
}

ParticleSystem::generateRandomClusters(float xmin, float xmax, float ymin, float ymax, float R){

}
