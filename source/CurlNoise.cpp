//
// CurlNoise.cpp
// Copyright (c) 2016-2018
// author: Douglas Creel
//

#include "CurlNoise.h"
#include "Utilities.h"

using namespace Utilities;

//
// public methods
//

CurlNoise::CurlNoise() : m_flow(0), m_deltaX(1e-4), m_lacunarity(1.0), m_gain(0.5)
{
	//m_noise.GenerateNoiseTile(256, 0);
}

CurlNoise::~CurlNoise()
{
}

void CurlNoise::setDeltaX(float dx)
{
	m_deltaX = dx;
}

// compute curl of potential field
vector3 CurlNoise::getCurlVelocity(vector3 pos, vector3 center, float flow)
{
	vector3 curl;
	pos.x -= center.x;
	pos.y -= center.y;
	pos.z -= center.z;
	pos.x *= m_frequency.x;
	pos.y *= m_frequency.y;
	pos.z *= m_frequency.z;
	m_flow = flow,

	curl.x = ((potential(pos.x, pos.y + m_deltaX, pos.z).z - potential(pos.x, pos.y - m_deltaX, pos.z).z) -
		(potential(pos.x, pos.y, pos.z + m_deltaX).y - potential(pos.x, pos.y, pos.z - m_deltaX).y)) / (2.0f * m_deltaX);
	curl.y = ((potential(pos.x, pos.y, pos.z + m_deltaX).x - potential(pos.x, pos.y, pos.z - m_deltaX).x) -
		(potential(pos.x + m_deltaX, pos.y, pos.z).z - potential(pos.x - m_deltaX, pos.y, pos.z).z)) / (2.0f * m_deltaX);
	curl.z = ((potential(pos.x + m_deltaX, pos.y, pos.z).y - potential(pos.x - m_deltaX, pos.y, pos.z).y) -
		(potential(pos.x, pos.y + m_deltaX, pos.z).x - potential(pos.x, pos.y - m_deltaX, pos.z).x)) / (2.0f * m_deltaX);

	return curl;
}

void CurlNoise::setFrequency(vector3 freq)
{
	m_frequency = freq;
}

//
// private methods
//

// define potential field
vector3 CurlNoise::potential(float x, float y, float z)
{
	vector3 retval;

	// use offsets into same noise to get 3 components of vector field
	float pnt[3] = {x, y, z};
	float pnt_offset[3] = { 0.f, 10.f, -10.f };

	pnt[0] += pnt_offset[0]; pnt[1] += pnt_offset[0]; pnt[2] += pnt_offset[0];

	retval.x = m_noise.turbulence(pnt, m_lacunarity, m_gain, m_octaves, m_flow);

	pnt[0] += pnt_offset[1]; pnt[1] += pnt_offset[1]; pnt[2] += pnt_offset[1];

	retval.y = m_noise.turbulence(pnt, m_lacunarity, m_gain, m_octaves, m_flow);

	//pnt[0] += -pnt_offset[1] + pnt_offset[2]; pnt[1] += -pnt_offset[1] + pnt_offset[2]; pnt[2] += -pnt_offset[1] + pnt_offset[2];
	pnt[0] += pnt_offset[2]; pnt[1] += pnt_offset[2]; pnt[2] += pnt_offset[2];

	retval.z = m_noise.turbulence(pnt, m_lacunarity, m_gain, m_octaves, m_flow);

	return retval;
}

void CurlNoise::setTurbulenceSettings(float lacunarity, float gain, int octaves)
{
	m_lacunarity = lacunarity;
	m_gain = gain;
	m_octaves = octaves;
}