#include "../headers/pointcloud.h"
using namespace std;
#include "../utils/simple_math.h"
#include "../utils/simple_timer.h"

#include "../headers/hashtable.h"
#include "../headers/pointcloud.h"

#include <unordered_map>



Grid par;


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










long long int mulmod(long long int a, long long int b, long long int m) {//It returns true if number is prime otherwise false {
   long long int x = 0, y = a % m;
   while (b > 0) {
      if (b % 2 == 1) {
         x = (x + y) % m;
      }
      y = (y * 2) % m;
      b /= 2;
   }
   return x % m;
}

long long int modulo(long long int base, long long int e, long long int m) {
   long long int x = 1;
   long long int y = base;
   while (e > 0) {
      if (e % 2 == 1)
         x = (x * y) % m;
         y = (y * y) % m;
         e = e / 2;
   }
   return x % m;
}

bool Miller(long long int p, int iteration) {
   if (p < 2) {
      return false;
   }
   if (p != 2 && p % 2==0) {
      return false;
   }
   long long int s = p - 1;
   while (s % 2 == 0) {
      s /= 2;
   }
   for (int i = 0; i < iteration; i++) {
      long long int a = rand() % (p - 1) + 1, temp = s;
      long long int mod = modulo(a, temp, p);
      while (temp != p - 1 && mod != 1 && mod != p - 1) {
         mod = mulmod(mod, mod, p);
         temp *= 2;
      }
      if (mod != p - 1 && temp % 2 == 0) {
         return false;
      }
   }
   return true;
}

long long int smallest_prime_atleast(long long int m){
	long long int n = m;
	while(!Miller(n, 10)) ++n;
	return n;
}


#include "fofcpp.cppx"

void PointCloud::group_grid_fof(float Rg){

//	SimpleTimer T; T.reset(); T.start();
//	calcGridParams(Rg/sqrt(3), par);
//	T.stop(); T.printTime("calcGrid");

	
	vector<double> pos_vec_d;
	pos_vec_d.assign(points.begin(), points.end());
	
	double * pos = pos_vec_d.data();


	int num_pos = nverts;
	
	float rcut = Rg;
	
	double minmax[6];	
	get_min_max(pos, num_pos, minmax);
	
	float Lx = minmax[3]-minmax[0], Ly=minmax[4]-minmax[1], Lz=minmax[5]-minmax[2];
	cout << "minmax: " << minmax[3] << "-" << minmax[0] << "," << minmax[4] << "-" << minmax[1] << "," <<  minmax[5] << "-" << minmax[2] << endl;
	cout << "dx dy dz = " << Lx << " " << Ly << " " << Lz << endl;
	float max_dim = fmax(fmax(Lx, Ly), Lz);
	
	
	float cell_width = rcut / sqrt(3);
	float inv_cell_width = 1.0/cell_width;

    int n_min = (int)(ceil(max_dim*inv_cell_width))+2; // two extra cells so never wrap

	int N = smallest_prime_atleast(n_min); // Make sure a prime
	int M = smallest_prime_atleast(N*N);

	printf("n_min = %d, N = %d, N2 = %d, M = %d\n", n_min, N, N*N, M);

	
	vector<int64_t> cells(num_pos);
	find_lattice(pos, num_pos, inv_cell_width, N, M, cells.data());
	
//	for (int i=0; i<num_pos; ++i) cout << cells[i] << " ";
//	cout << endl;

	vector <int64_t> sort_idx(num_pos); // = {0,1,2,3,4}; // to be found by index_sorting
	for (int i=0; i<num_pos; ++i) sort_idx[i] = i;

	SimpleTimer T; 
	T.reset(); T.start();
	sort(sort_idx.begin(), sort_idx.end(), [&cells](unsigned int i, unsigned int j){return cells[i] < cells[j];}); 
	T.stop(); T.printTime("sorting");

//	for (int i=0; i<num_pos; ++i) cout << sort_idx[i] << " ";
//	cout << endl;

	
	T.reset(); T.start();
	vector <int32_t> domains(num_pos);
	fof_link_cells(num_pos, N, M, rcut, 
		   pos, cells.data(), 
		   sort_idx.data(), domains.data());
		   
//	for (int i=0; i<num_pos; ++i) cout << domains[i] << " ";
//	cout << endl;
	T.stop(); T.printTime("fof");
	
	
	group_ids.assign(domains.begin(), domains.end());
	
		

}
	




#include "fof64.cppx"

void PointCloud::group_grid_fof64(float Rg){

//	SimpleTimer T; T.reset(); T.start();
//	calcGridParams(Rg/sqrt(3), par);
//	T.stop(); T.printTime("calcGrid");

	
	vector<double> pos_vec_d;
	pos_vec_d.assign(points.begin(), points.end());
	
	double * pos = pos_vec_d.data();


	int num_pos = nverts;
	
	float rcut = Rg;
	
	double minmax[6];	
	get_min_max(pos, num_pos, minmax);
	
	float Lx = minmax[3]-minmax[0], Ly=minmax[4]-minmax[1], Lz=minmax[5]-minmax[2];
	cout << "minmax: " << minmax[3] << "-" << minmax[0] << "," << minmax[4] << "-" << minmax[1] << "," <<  minmax[5] << "-" << minmax[2] << endl;
	cout << "dx dy dz = " << Lx << " " << Ly << " " << Lz << endl;
	float max_dim = fmax(fmax(Lx, Ly), Lz);
	

	float cell_width = rcut / sqrt(3);
	float inv_cell_width = 1.0/cell_width;


    float inv_block_width = 0.25 * inv_cell_width;
    
    int n_min = (int)(ceil(max_dim*inv_block_width))+2; // two extra cells so never wrap

	int N = smallest_prime_atleast(n_min); // Make sure a prime
	int M = smallest_prime_atleast(N*N);

	printf("n_min = %d, N = %d, N2 = %d, M = %d\n", n_min, N, N*N, M);

	
	vector<int64_t> blockcells(num_pos);
	blocks_cells(minmax[0], minmax[1], minmax[2], 
		  pos, num_pos, 
		  inv_cell_width, N, M, 
		  blockcells.data());

//	for (int i=0; i<num_pos; ++i) cout << blockcells[i] << " ";
//	cout << endl;

	vector <int64_t> sort_idx(num_pos); // = {0,1,2,3,4}; // to be found by index_sorting
	for (int i=0; i<num_pos; ++i) sort_idx[i] = i;

	SimpleTimer T; 
	T.reset(); T.start();
	sort(sort_idx.begin(), sort_idx.end(), [&blockcells](unsigned int i, unsigned int j){return blockcells[i] < blockcells[j];}); 
	T.stop(); T.printTime("sorting");

//	for (int i=0; i<num_pos; ++i) cout << sort_idx[i] << " ";
//	cout << endl;


	T.reset(); T.start();
	vector <int32_t> domains(num_pos);
	fof64(num_pos, N, M, num_pos, rcut, 
		  pos, blockcells.data(), 
		  sort_idx.data(), NULL,
		  domains.data(), 0.6);
	  
	  		   
//	for (int i=0; i<num_pos; ++i) cout << domains[i] << " ";
//	cout << endl;
	T.stop(); T.printTime("fof");
	
	
	group_ids.assign(domains.begin(), domains.end());
	
		

}
	











