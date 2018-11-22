#ifndef PARTICLES_H
#define PARTICLES_H

#include <string>
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
//#include <cuda_runtime.h>
using namespace std;

class ParticleSystem{
	public:
	string name;

	int N;
//	vector <Particle> pvec;

	float * pos, pos_dev;

	public:

//	ParticleSystem(int N);
//	~ParticleSystem();

	void generateRandomClusters(float xmin, float xmax, float ymin, float ymax, float R, int nc);
	
	void updateGroupIndices_grid();

};



#endif


