//
// SimplexRDNoise.h 
//
// code attribution: Stefan Gustavson https://github.com/stegu/perlin-noise (from original Simplex noise by Ken Perlin)
//

#ifndef SIMPLEXRDNOISE_H
#define SIMPLEXRDNOISE_H

#include <math.h>
#include <algorithm> 
#include "Utilities.h"

using std::floor;


class SimplexRDNoise
{
public:
	SimplexRDNoise::SimplexRDNoise();

	SimplexRDNoise::~SimplexRDNoise();

	float SimplexRDNoise::SimplexRDNoise2D(float x, float y, float angle, float *dnoise_dx, float *dnoise_dy);

	float SimplexRDNoise::SimplexRDNoise3D(float x, float y, float z, float angle, float *dnoise_dx, float *dnoise_dy, float *dnoise_dz);

	float SimplexRDNoise::turbulence(float pnt[3], float lacunarity, float gain, int octaves, float flow);

private:
	void gradrot2(int hash, float sin_t, float cos_t, float *gx, float *gy);

	void gradrot3(int hash, float sin_t, float cos_t, float *gx, float *gy, float *gz);

	float graddotp2(float gx, float gy, float x, float y);

	float graddotp3(float gx, float gy, float gz, float x, float y, float z);
};

#endif
