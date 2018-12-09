#ifndef HASHTABLE_CPU_GPU
#define HASHTABLE_CPU_GPU
#include <iostream>
using namespace std;

#include "../utils/simple_math.h"

struct KeyHasher{
	inline int operator()(const int3& key) const {
		return (key.z*377771 + key.y*677 + key.x);	// 2.15, 64
	}
};


// HASH TABLE IMPLEMENTATION OF GRID

class HashNode{
	public:
	int3 key;
	int2 value;
	int count;
	
	inline DEVICE_NAME HashNode(){
		key.x = key.y = key.z = -1;
		value.x = value.y = 0;
		count = 0;
	}
};


// HASH TABLE IMPLEMENTATION OF GRID


inline DEVICE_NAME int hash3D(int3 key, int length){
	return (key.z*377771 + key.y*677 + key.x) % length;	// 2.15, 64
//	return (key.z*377771 + key.y*133337 + key.x) % length; // 1.7, 25
//	return (key.z*497771 + key.y*133337 + key.x) % length; // 30, 173
//	return ((long long int)key.z*497771 + key.y*787 + key.x) % length; // 1.5, 28
}


inline DEVICE_NAME int hash_insert(int3 key, int2 value, HashNode* ht, int length, int* final_id=NULL	){
	int id = hash3D(key, length);
	int hash = id;
	int count = 1;
	while (ht[id].key != HashNode().key && ht[id].key != key){ 
		id = (id+11)%length;
		++count;
	}
	ht[id].value = value;
	ht[id].key = key;
	ht[hash].count = count; //max(counts[id], count);
	if (final_id != NULL) *final_id = id;
	return count;

}


inline DEVICE_NAME int hash_find(int3 key, HashNode* ht, int length, int* attempts=NULL){
	int id = hash3D(key, length);
	int hash = id;
	int count = 1; 
	if (attempts != NULL) *attempts = count;
	while (ht[id].key != key){ //(keys[id].x != key.x || keys[id].y != key.y || keys[id].z != key.z){
		id = (id+11)%length;
		++count;
		if (attempts != NULL) *attempts = count;
		if (count > ht[hash].count) return -1;
	}
	return id;
}


#endif

