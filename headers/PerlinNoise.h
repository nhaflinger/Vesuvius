//
// PerlinNoise.h 
//
// code attribution: Stefan Gustavson (from original Perlin noise by Ken Perlin)
//

#ifndef PERLINNOISE_H
#define PERLINNOISE_H

#include <math.h>
#include <algorithm> 
#include "Utilities.h"

using std::floor;


class PerlinNoise
{
public:
	PerlinNoise();

	~PerlinNoise();

	float PerlinNoise1D(float x);

	float PerlinNoise2D(float x, float y);

	float PerlinNoise3D(float x, float y, float z);

	float PerlinNoise4D(float x, float y, float z, float w);

	float turbulence(float pnt[3], float lacunarity, float gain, int octaves, float flow);

private:
	float fade(float t);

	float grad1D(int hash, float x);

	float grad2D(int hash, float x, float y);

	float grad3D(int hash, float x, float y, float z);

	float grad4D(int hash, float x, float y, float z, float w);
};


#endif
