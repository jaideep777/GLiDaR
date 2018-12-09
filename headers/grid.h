#ifndef GRID_CPU_GPU
#define GRID_CPU_GPU

#include "../utils/simple_math.h"

class Grid{
	public:
	float3 cellSize;
	int3 gridSize;
	float3 origin;
	
	inline DEVICE_NAME int3 pos2cell(float3 pos){
		int3 cell;
		cell.x = (pos.x - origin.x)/cellSize.x;
		cell.y = (pos.y - origin.y)/cellSize.y;
		cell.z = (pos.z - origin.z)/cellSize.z;
		return cell;
	}

	inline DEVICE_NAME int3 index2cell(int cellID){
		int3 cell;
		cell.z = int( cellID / gridSize.x / gridSize.y);
		cell.y = int((cellID - cell.z*gridSize.x*gridSize.y) / gridSize.x);
		cell.x = int((cellID - cell.z*gridSize.x*gridSize.y - cell.y*gridSize.x));
		return cell;
	}

	inline DEVICE_NAME int cell2index(int3 cell){
		//		   iz * nx * ny             +    iy * nx        + ix
		return cell.z*gridSize.x*gridSize.y + cell.y*gridSize.x + cell.x;
		// can use z-order curve as hash here rather than 1D cell-index
	}

};


#endif
