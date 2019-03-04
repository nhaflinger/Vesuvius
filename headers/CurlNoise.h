//
// CurlNoise.h 
// Copyright (c) 2016-2018
// author: Douglas Creel
//

#ifndef CURLNOISE_H
#define CURLNOISE_H

#include "SimplexNoise.h"


class CurlNoise
{
public:
	CurlNoise();

	~CurlNoise();

	vector3 getCurlVelocity(vector3 pos, vector3 center, float flow);

	void setDeltaX(float dx);
	
	void setTurbulenceSettings(float lacunarity, float gain, int octaves);

	void setFrequency(vector3 freq);

private:
	vector3 potential(float x, float y, float z);

	float m_deltaX;

	float m_flow;

	SimplexNoise m_noise;

	float m_lacunarity = 2.0f;
	
	float m_gain = 0.5f;
	
	int m_octaves = 1;

	vector3 m_frequency = { 1.0f, 1.0f, 1.0f };
};


#endif
