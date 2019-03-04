//
// WaveletNoise.h 
// Copyright (c) 2016-2018
// author: Douglas Creel
//
// code attribution: "Wavelet Noise" by Robert L. Cook & Tony DeRose
//

#ifndef WAVELETNOISE_H
#define WAVELETNOISE_H

#include "Vesuvius.h"
#include <math.h>
#include <random>
#include <Eigen/Dense>


#define ARAD 16 



class WaveletNoise
{
public:
	WaveletNoise();

	WaveletNoise(int noiseTileSize);

	~WaveletNoise();

	float WNoise(float p[3]);

	float WProjectedNoise(float p[3], float normal[3]);

	float WMultibandNoise(float p[3], float s, float *normal, int firstBand, int nbands, float *w);

	void GenerateNoiseTile(int n, int olap);

	float gaussianNoise(float mean, float stdDev);

	void setRandSeed(int seed) { m_randseed = seed; }

	float fBm(float pnt[3], float lacunarity, float gain, int octaves);

	float turbulence(float pnt[3], float lacunarity, float gain, int octaves, float flow);


private:
	void Downsample(float *from, float *to, int n, int stride);

	void Upsample(float *from, float *to, int n, int stride);

	int m_noiseTileSize = 128;

	std::vector<float> m_noiseTileData;

	int m_randseed = 0;
};

#endif

