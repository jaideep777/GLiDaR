#include "../headers/hashtable.h"
using namespace std;

#ifdef __CUDACC__

#define DEVICE_NAME __device__ __host__

#else

#define DEVICE_NAME 

bool int3::operator != (int3 v){
	return (x != v.x || y != v.y || z != v.z);
}

bool int3::operator==(const int3 &other) const{ 
	return (x == other.x && y== other.y && z == other.z);
}

#endif


bool DEVICE_NAME compare(int3 i, int3 j){
	return (i.x == j.x && i.y == j.y && i.z == j.z);
}

// HASH TABLE IMPLEMENTATION OF GRID

DEVICE_NAME HashNode::HashNode(){
	key.x = key.y = key.z = -1;
	value.x = value.y = 0;
	count = 0;
}



int DEVICE_NAME hash3D(int3 key, int length){
	return (key.z*377771 + key.y*677 + key.x) % length;	// 2.15, 64
//	return (key.z*377771 + key.y*133337 + key.x) % length; // 1.7, 25
//	return (key.z*497771 + key.y*133337 + key.x) % length; // 30, 173
//	return ((long long int)key.z*497771 + key.y*787 + key.x) % length; // 1.5, 28
}


int DEVICE_NAME hash_insert(int3 key, int2 value, HashNode* ht, int length, int* final_id){
	int id = hash3D(key, length);
	int hash = id;
	int count = 1;
	while (!compare(ht[id].key, HashNode().key) && !compare(ht[id].key, key)){ 
		id = (id+11)%length;
		++count;
	}
	ht[id].value = value;
	ht[id].key = key;
	ht[hash].count = count; //max(counts[id], count);
	if (final_id != NULL) *final_id = id;
	return count;

}


int DEVICE_NAME hash_find(int3 key, HashNode* ht, int length, int* attempts){
	int id = hash3D(key, length);
	int hash = id;
	int count = 1; 
	if (attempts != NULL) *attempts = count;
	while (!compare(ht[id].key, key)){ //(keys[id].x != key.x || keys[id].y != key.y || keys[id].z != key.z){
		id = (id+11)%length;
		++count;
		if (attempts != NULL) *attempts = count;
		if (count > ht[hash].count) return -1;
	}
	return id;
}



