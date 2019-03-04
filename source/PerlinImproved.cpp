//
// PerlinImproved.cpp
// Copyright (c) 2016 Pixel Grammar, LLC
// author: Douglas Creel
//

#include "PerlinImproved.h"
#include <iostream>

using namespace Utilities;


//
// public methods
//

PerlinImproved::PerlinImproved()
{
	for (int x = 0; x<512; x++) 
	{
		perm[x] = permutation[x % 256];
	}
}

PerlinImproved::~PerlinImproved()
{
}

double PerlinImproved::ImprovedNoise(double x, double y, double z, double flow)
{
	int px = (int)floor(x) & 255;
	int py = (int)floor(y) & 255;
	int pz = (int)floor(z) & 255;

	double nx = x - floor(x);                               
	double ny = y - floor(y);
	double nz = z - floor(z);

	double u = fade(nx);
	double v = fade(ny);
	double w = fade(nz);

	int A = perm[px] + py;
	int AA = perm[A] + pz;
	int AB = perm[A + 1] + pz;
	int B = perm[px + 1] + py;
	int BA = perm[B] + pz;
	int BB = perm[B + 1] + pz;

	double noiseval = lerp(w, lerp(v, lerp(u, grad(perm[AA], nx, ny, nz),
		                                      grad(perm[BA], nx - 1, ny, nz)),
		                              lerp(u, grad(perm[AB], nx, ny - 1, nz),
		                                      grad(perm[BB], nx - 1, ny - 1, nz))),
		                      lerp(v, lerp(u, grad(perm[AA + 1], nx, ny, nz - 1),
		                                      grad(perm[BA + 1], nx - 1, ny, nz - 1)),
		                              lerp(u, grad(perm[AB + 1], nx, ny - 1, nz - 1),
	                                          grad(perm[BB + 1], nx - 1, ny - 1, nz - 1))));

	return noiseval;
}

double PerlinImproved::turbulence(double pnt[3], double lacunarity, double gain, int octaves, double flow)
{
	double amp = 1.0;
	double sum = 0.0f;
	double pp[3]; pp[0] = pnt[0]; pp[1] = pnt[1]; pp[2] = pnt[2];

	for (int i = 0; i < octaves; i++)
	{
		sum += amp * ImprovedNoise(pp[0], pp[1], pp[2], flow);
		amp *= gain;
		pp[0] *= lacunarity; pp[1] *= lacunarity; pp[2] *= lacunarity;
	}

	return sum;
}


void PerlinImproved::setFlowPulse(double pulse)
{
	m_flowPulse = pulse;
}

//
// private methods
//

double PerlinImproved::fade(double t)
{
	return t * t * t * (t * (t * 6 - 15) + 10);
}
/*
double PerlinImproved::grad(int hash, double x, double y, double z)
{
	int h = hash & 15;
	double u = h < 8 ? x : y;
	double v = h<4 ? y : h == 12 || h == 14 ? x : z;
	return ((h & 1) == 0 ? u : -u) + ((h & 2) == 0 ? v : -v);
}
*/
// easier to understand version of hash function
double PerlinImproved::grad(int hash, double x, double y, double z)
{
	switch (hash & 0xF)
	{
	case 0x0: return  x + y;
	case 0x1: return -x + y;
	case 0x2: return  x - y;
	case 0x3: return -x - y;
	case 0x4: return  x + z;
	case 0x5: return -x + z;
	case 0x6: return  x - z;
	case 0x7: return -x - z;
	case 0x8: return  y + z;
	case 0x9: return -y + z;
	case 0xA: return  y - z;
	case 0xB: return -y - z;
	case 0xC: return  y + x;
	case 0xD: return -y + z;
	case 0xE: return  y - x;
	case 0xF: return -y - z;
	default: return 0; 
	}
}