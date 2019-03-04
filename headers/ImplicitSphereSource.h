//
// ImplicitSphereSource.h 
// Copyright (c) 2016-2018
// author: Douglas Creel
//

#ifndef IMPLICITSPHERESOURCE_H
#define IMPLICITSPHERESOURCE_H


#include "Vesuvius.h" 
#include "Grid3D.h"
#include "FluidSource.h"
#include "PerlinNoise.h"
#include "SimplexNoise.h"
#include "WaveletNoise.h"


class ImplicitSphereSource : public FluidSource
{
public:
    ImplicitSphereSource(bool activate, UPDATEBEHAVIOR behavior, float featherAmount, bool emptyInterior, float dx, float dt, vector3 center, float radius, float density, Grid3D<float>& smokeGrid, float densityScale, float temperature, Grid3D<float>& temperatureGrid, float temperatureScale, float fuel, Grid3D<float>& fuelGrid, float fuelScale, bool enableNoise, NoiseParams* noiseParam, SimplexNoise noise, vector3 velocity, Grid3D<float>& veluGrid, Grid3D<float>& velvGrid, Grid3D<float>& velwGrid);

	~ImplicitSphereSource();

	virtual void update(Grid3D<float>& smokeGrid, Grid3D<float>& temperatureGrid, Grid3D<float>& fuelGrid, Grid3D<float>& veluGrid, Grid3D<float>& velvGrid, Grid3D<float>& velwGrid);

	void setFeatherAmount(float featherAmount);

	void enableNoise(bool enable);

	void setNoiseAmpScale(float scale);

	void setNoiseFreqScale(vector3 scale);

	void setNoisePhaseOffset(float offset);

	void setCenter(vector3 center);

	void setRadius(float radius);

	void setDensity(float density);

	void setTemperature(float temperature);

	void setFuel(float fuel);

	void setDensityScale(float scale);

	void setTemperatureScale(float scale);

	void setFuelScale(float scale);

	void setTimeStep(float step);

	void setVoxelSize(float voxelsize);

	void setVelocity(vector3 velocity);

private:
	void setVeluGrid(Grid3D<float>& veluGrid);

	void setVelvGrid(Grid3D<float>& velvGrid);

	void setVelwGrid(Grid3D<float>& velwGrid);

	float m_dx;

	float m_dt;

	vector3 m_center; 

	float m_radius;

	float m_density;

	float m_temperature;

	float m_fuel;

	float m_densityScale = 1.0f;

	float m_temperatureScale = 1.0f;

	float m_fuelScale = 1.0f;
	
	float m_featherAmount = 0.2f;

	bool m_enableNoise = false;

	NoiseParams* m_noiseParamStruct;

	bool m_emptyInterior = true;

	SimplexNoise m_noise;

	vector3 m_velocity = { 0, 0, 0 };
};
#endif
