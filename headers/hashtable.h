#ifndef HASHTABLE_CPU_GPU
#define HASHTABLE_CPU_GPU
#include <iostream>
using namespace std;

#ifdef __CUDACC__

#define DEVICE_NAME __device__ __host__

#else

#define DEVICE_NAME 

class int2{
	public:
	int x, y;
};

class int3{
	public:
	int x;
	int y;
	int z;
	
	bool operator != (int3 v);
	bool operator==(const int3 &other) const;

};

#endif


bool DEVICE_NAME compare(int3 i, int3 j);

struct KeyHasher{
	int operator()(const int3& key) const {
		return (key.z*377771 + key.y*677 + key.x);	// 2.15, 64
	}
};


// HASH TABLE IMPLEMENTATION OF GRID

class HashNode{
	public:
	int3 key;
	int2 value;
	int count;
	
	HashNode();
};


int DEVICE_NAME hash3D(int3 key, int length);

int DEVICE_NAME hash_insert(int3 key, int2 value, HashNode* ht, int length, int* final_id = NULL);

int DEVICE_NAME hash_find(int3 key, HashNode* ht, int length, int* attempts = NULL);


#endif

