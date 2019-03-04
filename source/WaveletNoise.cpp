//
// WaveletNoise.cpp
// Copyright (c) 2016-2018
//
// Implementation by Douglas Creel
// 
// Taken from "Wavelet Noise" by Robert L. Cook and Tony DeRose
//

#include "WaveletNoise.h"
#include <iostream>

using namespace Utilities;

int Mod(int x, int n) { int m = x%n; return (m<0) ? m + n : m; }


//
// public methods
//

WaveletNoise::WaveletNoise()
{
}

WaveletNoise::WaveletNoise(int noiseTileSize)
{
	GenerateNoiseTile(noiseTileSize, 0);
}

WaveletNoise::~WaveletNoise()
{
}

float WaveletNoise::WNoise(float p[3])  /* Non-projected 3D noise */
{
	int i, f[3], c[3], mid[3], n = m_noiseTileSize; /* f, c = filter, noise coeff indices */
	float w[3][3], t, result = 0;

	/* Evaluate quadratic B-spline basis functions */
	for (i = 0; i<3; i++) 
	{
		mid[i] = ceil(p[i] - 0.5); t = mid[i] - (p[i] - 0.5);
		w[i][0] = t*t / 2; w[i][2] = (1 - t)*(1 - t) / 2; w[i][1] = 1 - w[i][0] - w[i][2];
	}

	/* Evaluate noise by weighting noise coefficients by basis function values */
	for (f[2] = -1; f[2] <= 1; f[2]++) 
	{
		for (f[1] = -1; f[1] <= 1; f[1]++)
		{
			for (f[0] = -1; f[0] <= 1; f[0]++)
			{
				float weight = 1;

				for (i = 0; i < 3; i++)
				{
					c[i] = Mod(mid[i] + f[i], n);
					weight *= w[i][f[i] + 1];
				}
				result += weight * m_noiseTileData[c[2] * n*n + c[1] * n + c[0]];
			}
		}
	}
	
	return result;
}

float WaveletNoise::WProjectedNoise(float p[3], float normal[3]) /* 3D noise projected onto 2D */
{
	int i, c[3], min[3], max[3], n = m_noiseTileSize; /* c = noise coeff location */
	float support, result = 0;

	/* Bound the support of the basis functions for this projection direction */
	for (i = 0; i<3; i++)
	{
		support = 3 * abs(normal[i]) + 3 * sqrt((1 - normal[i] * normal[i]) / 2);
		min[i] = ceil(p[i] - (3 * abs(normal[i]) + 3 * sqrt((1 - normal[i] * normal[i]) / 2)));
		max[i] = floor(p[i] + (3 * abs(normal[i]) + 3 * sqrt((1 - normal[i] * normal[i]) / 2)));
	}

	/* Loop over the noise coefficients within the bound. */
	for (c[2] = min[2]; c[2] <= max[2]; c[2]++)
	{
		for (c[1] = min[1]; c[1] <= max[1]; c[1]++)
		{
			for (c[0] = min[0]; c[0] <= max[0]; c[0]++)
			{
				float t, t1, t2, t3, dot = 0, weight = 1;
				/* Dot the normal with the vector from c to p */
				for (i = 0; i<3; i++) { dot += normal[i] * (p[i] - c[i]); }
				/* Evaluate the basis function at c moved halfway to p along the normal. */
				for (i = 0; i<3; i++)
				{
					t = (c[i] + normal[i] * dot / 2) - (p[i] - 1.5); t1 = t - 1; t2 = 2 - t; t3 = 3 - t;
					weight *= (t <= 0 || t >= 3) ? 0 : (t<1) ? t*t / 2 : (t<2) ? 1 - (t1*t1 + t2*t2) / 2 : t3*t3 / 2;
				}
				/* Evaluate noise by weighting noise coefficients by basis function values. */
				result += weight * m_noiseTileData[Mod(c[2], n)*n*n + Mod(c[1], n)*n + Mod(c[0], n)];
			}
		}
	}

	return result;
}

float WaveletNoise::WMultibandNoise(float p[3], float s, float *normal, int firstBand, int nbands, float *w)
{
	float q[3], result = 0, variance = 0; int i, b;

	for (b = 0; b<nbands && s + firstBand + b<0; b++)
	{
		for (i = 0; i <= 2; i++)
		{
			q[i] = 2 * p[i] * pow(2, firstBand + b);
		}
		result += (normal) ? w[b] * WProjectedNoise(q, normal) : w[b] * WNoise(q);
	}

	for (b = 0; b<nbands; b++)
	{
		variance += w[b] * w[b];
	}

	/* Adjust the noise so it has a variance of 1. */
	if (variance) result /= sqrt(variance * ((normal) ? 0.296 : 0.210));

	return result;
}

void WaveletNoise::GenerateNoiseTile(int n, int olap)
{
	if (n % 2) n++; /* tile size must be even */
	int ix, iy, iz, i, sz = n*n*n * sizeof(float);
	std::vector<float> temp1; std::vector<float> temp2;
	temp1.reserve(sz); temp2.reserve(sz);
	m_noiseTileData.reserve(sz);

	/* Step 1. Fill the tile with random numbers in the range -1 to 1. */
	std::vector<float>::iterator it;
	srand(1001 + m_randseed);
	for (int i = 0; i < sz; i++) {
		m_noiseTileData.push_back(gaussianNoise(0, 1));
	}

	/* Steps 2 and 3. Downsample and upsample the tile */
	for (iy = 0; iy<n; iy++) for (iz = 0; iz<n; iz++) /* each x row */
	{
		i = iy*n + iz*n*n;
		Downsample(&m_noiseTileData[i], &temp1[i], n, 1);
		Upsample(&temp1[i], &temp2[i], n, 1);
	}

	for (ix = 0; ix<n; ix++) for (iz = 0; iz<n; iz++) /* each y row */
	{
		i = ix + iz*n*n;
		Downsample(&temp2[i], &temp1[i], n, n);
		Upsample(&temp1[i], &temp2[i], n, n);
	}

	for (ix = 0; ix<n; ix++) for (iy = 0; iy<n; iy++) /* each z row */
	{
		i = ix + iy*n;
		Downsample(&temp2[i], &temp1[i], n, n*n);
		Upsample(&temp1[i], &temp2[i], n, n*n);
	}

	/* Step 4. Subtract out the coarse-scale contribution */
	for (i = 0; i<n*n*n; i++) { m_noiseTileData[i] -= temp2[i]; }

	/* Avoid even/odd variance difference by adding odd-offset version of noise to itself.*/
	int offset = n / 2; if (offset % 2 == 0) offset++;
	for (i = 0, ix = 0; ix<n; ix++) for (iy = 0; iy<n; iy++) for (iz = 0; iz<n; iz++)
		temp1[i++] = m_noiseTileData[Mod(ix + offset, n) + Mod(iy + offset, n)*n + Mod(iz + offset, n)*n*n];

	for (i = 0; i<n*n*n; i++)
	{
		m_noiseTileData[i] += temp1[i];
	}
 
	m_noiseTileSize = n;
}

float WaveletNoise::gaussianNoise(float mean, float stdDev)
{
	float u, v, s;
	do
	{
		u = ((float)rand() / (float)RAND_MAX) * 2.0 - 1.0;
		v = ((float)rand() / (float)RAND_MAX) * 2.0 - 1.0;
		s = u * u + v * v;
	} while (s >= 1 || s == 0);

	float mul = sqrt(-2.0 * log(s) / s);

	return mean + stdDev * u * mul;
}

float WaveletNoise::fBm(float pnt[3], float lacunarity, float gain, int octaves)
{
	float amp = 1.0;
	float sum = 0.0f;
	float pp[3]; pp[0] = pnt[0]; pp[1] = pnt[1]; pp[2] = pnt[2];

	for (int i = 0; i < octaves; i++)
	{
		sum += amp * 0.5 * (WNoise(pp) + 1.0);
		amp *= gain;
		pp[0] *= lacunarity; pp[1] *= lacunarity; pp[2] *= lacunarity;
	}

	return sum;
}

float WaveletNoise::turbulence(float pnt[3], float lacunarity, float gain, int octaves, float flow)
{
	float amp = 1.0;
	float sum = 0.0f;
	float pp[3] = { pnt[0], pnt[1], pnt[2] };

	for (int i = 0; i < octaves; i++)
	{
		sum += amp * WNoise(pp);
		amp *= gain;
		pp[0] *= lacunarity; pp[1] *= lacunarity; pp[2] *= lacunarity;
	}

	return sum;
}

//
// private methods
//

void WaveletNoise::Downsample(float *from, float *to, int n, int stride)
{
	float *a, aCoeffs[2 * ARAD] = { 0.000334,-0.001528, 0.000410, 0.003545,-0.000938,-0.008233, 0.002172, 0.019120, -0.005040,-0.044412, 0.011655, 0.103311,-0.025936,-0.243780, 0.033979, 0.655340, 0.655340, 0.033979,-0.243780,-0.025936, 0.103311, 0.011655,-0.044412,-0.005040, 0.019120, 0.002172,-0.008233,-0.000938, 0.003546, 0.000410,-0.001528, 0.000334 };
	a = &aCoeffs[ARAD];

	for (int i = 0; i < n / 2; i++)
	{
		to[i*stride] = 0; 
		for (int k = 2 * i - ARAD; k <= 2 * i + ARAD; k++)
			to[i*stride] += a[k - 2 * i] * from[Mod(k, n)*stride];
	}
}

void WaveletNoise::Upsample(float *from, float *to, int n, int stride)
{
	float *p, pCoeffs[4] = { 0.25, 0.75, 0.75, 0.25 };
	p = &pCoeffs[2];

	for (int i = 0; i<n; i++)
	{
		to[i*stride] = 0;
		for (int k = i / 2; k <= i / 2 + 1; k++)
			to[i*stride] += p[i - 2 * k] * from[Mod(k, n / 2)*stride];
	}
}
