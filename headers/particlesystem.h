#ifndef PARTICLES_H
#define PARTICLES_H

#include <string>
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <cuda_runtime.h>


class ParticleSystem{
	public:
	string name;

	int N;
//	vector <Particle> pvec;

	float3 * pos, pos_dev;

	public:

	void generateRandomClusters(float xmin, float xmax, float ymin, float ymax, float R);
	
	void updateGroupIndices_grid();

};



#endif


