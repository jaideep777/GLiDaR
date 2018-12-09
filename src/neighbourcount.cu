#include <cuda_runtime.h>
#include <algorithm>
#include "../headers/hashtable.h"
#include "../headers/pointcloud.h"
using namespace std;

#include "../utils/simple_timer.h"
#include "../utils/cuda_vector_math.cuh"

float3 *points_dev;
int* pt_ids_dev;
HashNode* ht_dev;
int* nb_counts_dev;

__device__ void countNeighbours_hash_gpu(int3 c1, int3 c2, float Rg, int ht_size, float3*points_dev, int*pt_ids_dev, HashNode* ht_dev, int* nb_counts_dev){
	
	int count = 0;
	
//	int3 result; // (yes/no, point_id1, point_id2)
//	result.x = result.y = result.z = 0;
	
	int attempts1, attempts2;
	int id_c1 = hash_find(c1, ht_dev, ht_size, &attempts1);
	int id_c2 = hash_find(c2, ht_dev, ht_size, &attempts2);
	if (id_c1 == -1 || id_c2 == -1) return; // if either cell doesnt exist, cant merge cells
	
	for (int p1 = ht_dev[id_c1].value.x; p1 <= ht_dev[id_c1].value.y; ++p1){	
		for (int p2 = ht_dev[id_c2].value.x; p2 <= ht_dev[id_c2].value.y; ++p2){
			++count;
			
			int ip1 = pt_ids_dev[p1];
			int ip2 = pt_ids_dev[p2];
			
			float dx = points_dev[ip1].x - points_dev[ip2].x;
			float dy = points_dev[ip1].y - points_dev[ip2].y;
			float dz = points_dev[ip1].z - points_dev[ip2].z;
			
			float dist =  sqrt(dx*dx + dy*dy + dz*dz);

			if (dist < Rg){
				atomicAdd(&nb_counts_dev[ip1],1);
				atomicAdd(&nb_counts_dev[ip2],1);
			} 
		}
	}
//	cout << "attempts1 = " << attempts1 << " " << "attempts2 = " << attempts2 << endl;
	return; //result;
}


__global__ void countNeighbours_kernel(float Rd, float3*points_dev, int*pt_ids_dev, HashNode* ht_dev, int* nb_counts_dev, int ht_size, int npts, Grid par){
	
	int tid = blockIdx.x*blockDim.x + threadIdx.x;
	if (tid >= ht_size) return;

	int3 cell = ht_dev[tid].key;

	if (cell == HashNode().key) return;
	
//	int count=0;
	
	for (int ix = -2; ix <=2; ++ix){
		for (int iy = -2; iy <=2; ++iy){
			for (int iz = -2; iz <=2; ++iz){
				
				int3 cell_new;
				cell_new.x = clamp(cell.x + ix, 0, par.gridSize.x-1);
				cell_new.y = clamp(cell.y + iy, 0, par.gridSize.y-1);
				cell_new.z = clamp(cell.z + iz, 0, par.gridSize.z-1);

				if (par.cell2index(cell_new) <= par.cell2index(cell)){	
					countNeighbours_hash_gpu(cell, cell_new, Rd, ht_size, points_dev, pt_ids_dev, ht_dev, nb_counts_dev);
				}
			}
		}
	}
	

}


void PointCloud::countNeighbours_hash_gpu(float Rd){

	SimpleTimer T; T.reset(); T.start();
	Grid par;
	calcGridParams(Rd/sqrt(3), par);
	T.stop(); T.printTime("calcGrid");

	
	float3 * pos = (float3*)points.data();

	T.reset(); T.start();	
	vector <int> point_hashes(nverts);
	vector <int> point_ids(nverts);
	// get the cell ID for each particle
	for (int i=0; i<nverts; ++i) {
		int3 cell_id = par.pos2cell(pos[i]);
		point_hashes[i] = par.cell2index(cell_id);
		point_ids[i]    = i;
	}
	T.stop(); T.printTime("cell Ids");


	// sort particles by cell ID
	T.reset(); T.start();	
	sort(point_ids.begin(), point_ids.end(), [&point_hashes](int i, int j){return point_hashes[i] < point_hashes[j];}); 
	T.stop(); T.printTime("sort");


	int hashTable_size = 6000011;
	vector <HashNode> hashTable(hashTable_size);

	int n_attempts = 0;	
	int avg_attempts = 0;
	int max_attempts = 0;

	// Build Hashtable
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
		}
		++next;
	}
	cout << "Insertion Attempts: " << float(avg_attempts)/n_attempts << " " << max_attempts << endl;
	T.stop(); T.printTime("map cells");

	cudaMalloc(&ht_dev, hashTable_size*sizeof(HashNode));
	cudaMemcpy(ht_dev, hashTable.data(), hashTable_size*sizeof(HashNode), cudaMemcpyHostToDevice);
	
	cudaMalloc(&points_dev, 3*nverts*sizeof(float));
	cudaMemcpy(points_dev, points.data(), 3*nverts*sizeof(float), cudaMemcpyHostToDevice);

	cudaMalloc(&pt_ids_dev, nverts*sizeof(int));
	cudaMemcpy(pt_ids_dev, point_ids.data(), nverts*sizeof(int), cudaMemcpyHostToDevice);

	cudaMalloc(&nb_counts_dev, nverts*sizeof(int));


	neighbourCounts.resize(nverts,0);

//	particles_compared_per_cellpair = 0;
//	count_pcpc = 0;

//	// group grid cells
	int n_filled_cells = 0;
//	long long int pairs = 0;
	T.reset(); T.start();

	int nthreads = 256;
	int nblocks = (nverts-1)/nthreads+1;
	countNeighbours_kernel <<< nblocks, nthreads>>> (Rd, points_dev, pt_ids_dev, ht_dev, nb_counts_dev, hashTable_size, nverts, par);

	cudaMemcpy(neighbourCounts.data(), nb_counts_dev, nverts*sizeof(int), cudaMemcpyDeviceToHost);
	
//	for (int i=0; i<hashTable_size; ++i){
//		
//		int3 cell = hashTable[i].key;
//		if (compare(cell, HashNode().key)) continue;
//		++n_filled_cells;
//		
//		int count=0;
//		
//		for (int ix = -2; ix <=2; ++ix){
//			for (int iy = -2; iy <=2; ++iy){
//				for (int iz = -2; iz <=2; ++iz){
//					
//					int3 cell_new;
//					cell_new.x = clamp(cell.x + ix, 0, par.gridSize.x-1);
//					cell_new.y = clamp(cell.y + iy, 0, par.gridSize.y-1);
//					cell_new.z = clamp(cell.z + iz, 0, par.gridSize.z-1);
////					cout << "\tcell: " << cell_new.x << ", " << cell_new.y << ", " << cell_new.z << endl;

//					if (cellHash(cell_new) <= cellHash(cell)){	
////						int3 res = mergeCells(c1, c2, Rg, filled_cells, point_ids);
//						countNeighbours_hash(cell, cell_new, Rd, hashTable.data(), point_ids, hashTable_size, neighbourCounts);

//						++count;
//						++pairs;
//					}
////					if (pairs % 10000 == 0) cout << pairs << "pairs compared.\n"; 
//				}
//			}
//		}
//		
////		cells_compared_per_cell += count;
////		count_ccpc += 1;
//	}
//	cout << "pairs = " << pairs << endl;
	T.stop(); T.printTime("group cells");
	
	cout << "Summary: \n";
////	cout << "Particle pairs compared per mergeCell: " <<  float(particles_compared_per_cellpair)/count_pcpc << endl;
////	cout << "Cells compared per cell: " <<  float(cells_compared_per_cell)/count_ccpc << endl;
	cout << "No. of filled cells = " << n_filled_cells << endl;
	


}

