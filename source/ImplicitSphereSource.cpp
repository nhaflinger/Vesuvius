//
// ImplicitSphereSource.cpp
// Copyright (c) 2016-2018
// author: Douglas Creel
//


#include "ImplicitSphereSource.h"


using namespace Utilities;


ImplicitSphereSource::ImplicitSphereSource(bool activate, UPDATEBEHAVIOR behavior, float featherAmount, bool emptyInterior, float dx, float dt, vector3 center, float radius, float density, Grid3D<float>& smokeGrid, float densityScale, float temperature, Grid3D<float>& temperatureGrid, float temperatureScale, float fuel, Grid3D<float>& fuelGrid, float fuelScale, bool enableNoise, NoiseParams* noisestruct, SimplexNoise noise, vector3 velocity, Grid3D<float>& veluGrid, Grid3D<float>& velvGrid, Grid3D<float>& velwGrid)
{
	// create implicit sphere with a set smoke density and temperature
	m_active = activate;
	m_updateBehavior = behavior;
	m_center.x = center.x; m_center.y = center.y; m_center.z = center.z;
	m_radius = radius;
	m_density = density;
	m_temperature = temperature;
	m_fuel = fuel;
	m_dx = dx;
	m_dt = dt;
	m_enableNoise = enableNoise;
	m_emptyInterior = emptyInterior;
	m_featherAmount = featherAmount;
	m_noiseParamStruct = noisestruct;
	m_noise = noise;
	m_velocity = velocity;
	m_densityScale = densityScale;
	m_temperatureScale = temperatureScale;
	m_fuelScale = fuelScale;

#pragma omp parallel for    
	for (int i = 0; i < smokeGrid.getXSize(); i++)
	{
		for (int j = 0; j < smokeGrid.getYSize(); j++)
		{
			for (int k = 0; k < smokeGrid.getZSize(); k++)
			{
				vector3 dist;
				float x = i * dx; float y = j * dx; float z = k * dx;
				dist.x = x - center.x;
				dist.y = y - center.y;
				dist.z = z - center.z;
				float length = sqrt(dist.x*dist.x + dist.y*dist.y + dist.z*dist.z);
				float nudensity = density;
				float nutemperature = temperature;
				float nufuel = fuel;

				float tol = radius - length;
				if (m_emptyInterior)
					tol = abs(tol);
;
				if ((!m_emptyInterior && tol > 0.0) || (m_emptyInterior && tol < m_featherAmount))
				{
					if (m_enableNoise)
					{
						float flow = 0.0f;
						float pnt[3];
						pnt[0] = m_noiseParamStruct->noiseFreqScale.x * x;
						pnt[1] = (m_noiseParamStruct->noiseFreqScale.y * y + m_noiseParamStruct->noiseFlow);
						pnt[2] = m_noiseParamStruct->noiseFreqScale.z * z;
						float noise = m_noiseParamStruct->noiseAmpScale * m_noise.turbulence(pnt, m_noiseParamStruct->lacunarity, m_noiseParamStruct->gain, m_noiseParamStruct->octaves, m_noiseParamStruct->noiseFlow);
						// add noise
						if (m_noiseParamStruct->setMultAdd == MULTVALUE)
						{
							nudensity *= noise;
							nutemperature *= noise;
							nufuel *= noise;
						}
						else if (m_noiseParamStruct->setMultAdd == ADDVALUE)
						{
							nudensity += noise;
							nutemperature += noise;
							nufuel += noise;
						}
						else
						{
							nudensity = noise;
							nutemperature = noise;
							nufuel = noise;
						}
					}

					// feather input values
					if (!m_emptyInterior)
					{
						nudensity *= (1.0 - smoothstep(1.0 - m_featherAmount, 1.0, length / radius));
						nutemperature *= (1.0 - smoothstep(1.0 - m_featherAmount, 1.0, length / radius));
						nufuel *= (1.0 - smoothstep(1.0 - m_featherAmount, 1.0, length / radius));
					}
					else
					{
						nudensity *= (1.0 - smoothstep(0.0, 1.0, tol / m_featherAmount));
						nutemperature *= (1.0 - smoothstep(0.0, 1.0, tol / m_featherAmount));
						nufuel *= (1.0 - smoothstep(0.0, 1.0, tol / m_featherAmount));
					}

					nudensity *= m_densityScale;
					nutemperature *= m_temperatureScale;
					nufuel *= m_fuelScale;

					smokeGrid.set(nudensity, i, j, k);
					temperatureGrid.set(nutemperature, i, j, k);
					fuelGrid.set(nufuel, i, j, k);
				}
			}
		}
	}

	// update velocity grids
	setVeluGrid(veluGrid);
	setVelvGrid(velvGrid);
	setVelwGrid(velwGrid);
}

ImplicitSphereSource::~ImplicitSphereSource()
{
}

void ImplicitSphereSource::update(Grid3D<float>& smokeGrid, Grid3D<float>& temperatureGrid, Grid3D<float>& fuelGrid, Grid3D<float>& veluGrid, Grid3D<float>& velvGrid, Grid3D<float>& velwGrid)
{
	// update implicit sphere smoke density and temperature	

#pragma omp parallel for    
	for (int i = 0; i < smokeGrid.getXSize(); i++)
	{
		for (int j = 0; j < smokeGrid.getYSize(); j++)
		{
			for (int k = 0; k < smokeGrid.getZSize(); k++)
			{
				vector3 dist;
				float x = i * m_dx; float y = j * m_dx; float z = k * m_dx;
				dist.x = x - m_center.x;
				dist.y = y - m_center.y;
				dist.z = z - m_center.z;
				float length = sqrt(dist.x*dist.x + dist.y*dist.y + dist.z*dist.z);
				float nudensity = m_density;
				float nutemperature = m_temperature;
				float nufuel = m_fuel;

				float tol = m_radius - length;
				if (m_emptyInterior)
					tol = abs(tol);
				
				if ((!m_emptyInterior && tol > 0.0) || (m_emptyInterior && tol < m_featherAmount))
				{
					float smoke = smokeGrid(i, j, k);
					float temp = temperatureGrid(i, j, k);
					float fuel = fuelGrid(i, j, k);

					if (m_enableNoise)
					{
						float flow = 0.0f;
						float pnt[3];
						pnt[0] = m_noiseParamStruct->noiseFreqScale.x * x;
						pnt[1] = (m_noiseParamStruct->noiseFreqScale.y * y + m_noiseParamStruct->noiseFlow);
						pnt[2] = m_noiseParamStruct->noiseFreqScale.z * z;
						float noise = m_noiseParamStruct->noiseAmpScale * m_noise.turbulence(pnt, m_noiseParamStruct->lacunarity, m_noiseParamStruct->gain, m_noiseParamStruct->octaves, m_noiseParamStruct->noiseFlow);
						// add noise
						if (m_noiseParamStruct->setMultAdd == MULTVALUE)
						{
							nudensity *= noise;
							nutemperature *= noise;
							nufuel *= noise;
						}
						else if (m_noiseParamStruct->setMultAdd == ADDVALUE)
						{
							nudensity += noise;
							nutemperature += noise;
							nufuel += noise;
						}
						else
						{
							nudensity = noise;
							nutemperature = noise;
							nufuel = noise;
						}
					}

					// feather update values
					if (!m_emptyInterior)
					{
						nudensity *= (1.0 - smoothstep(1.0 - m_featherAmount * m_radius, 1.0, length / m_radius));
						nutemperature *= (1.0 - smoothstep(1.0 - m_featherAmount, 1.0, length / m_radius));
						nufuel *= (1.0 - smoothstep(1.0 - m_featherAmount, 1.0, length / m_radius));
					}
					else
					{
						nudensity *= (1.0 - smoothstep(0.0, 1.0, tol / m_featherAmount));
						nutemperature *= (1.0 - smoothstep(0.0, 1.0, tol / m_featherAmount));
						nufuel *= (1.0 - smoothstep(0.0, 1.0, tol / m_featherAmount));
					}

					nudensity *= m_densityScale;
					nutemperature *= m_temperatureScale;
					nufuel *= m_fuelScale;

					if (m_updateBehavior == COPY)
					{
						smokeGrid.set(nudensity, i, j, k);
						temperatureGrid.set(nutemperature, i, j, k);
						fuelGrid.set(nufuel, i, j, k);
					}
					else if (m_updateBehavior == ADD)
					{
						float nusmoke = smoke + m_dt * nudensity;
						smokeGrid.set(nusmoke, i, j, k);
						float nutemp = (temp + (1.0 - exp(-m_dt)) * (nutemperature - temp));
						temperatureGrid.set(nutemp, i, j, k);
						float nuf = fuel + m_dt * nufuel;
						fuelGrid.set(nuf, i, j, k);
					}
				}
			}
		}
	}

	// update velocity grids
	setVeluGrid(veluGrid);
	setVelvGrid(velvGrid);
	setVelwGrid(velwGrid);
}

void ImplicitSphereSource::setCenter(vector3 center)
{
	m_center = center;
}

void ImplicitSphereSource::setRadius(float radius)
{
	m_radius = radius;
}

void ImplicitSphereSource::setDensity(float density)
{
	m_density = density;
}

void ImplicitSphereSource::setTemperature(float temperature)
{
	m_temperature = temperature;
}

void ImplicitSphereSource::setDensityScale(float scale)
{
	m_densityScale = scale;
}

void ImplicitSphereSource::setTemperatureScale(float scale)
{
	m_temperatureScale = scale;
}

void ImplicitSphereSource::setFuelScale(float scale)
{
	m_fuelScale = scale;
}

void ImplicitSphereSource::setTimeStep(float step)
{
	m_dt = step;
}

void ImplicitSphereSource::setVoxelSize(float voxelsize)
{
	m_dx = voxelsize;
}

void ImplicitSphereSource::setFeatherAmount(float featherAmount)
{
	m_featherAmount = featherAmount;
}

void ImplicitSphereSource::enableNoise(bool enable)
{
	m_enableNoise = enable;
}

void ImplicitSphereSource::setNoiseAmpScale(float scale)
{
	m_noiseParamStruct->noiseAmpScale = scale;
}

void ImplicitSphereSource::setNoiseFreqScale(vector3 scale)
{
	m_noiseParamStruct->noiseFreqScale = scale;
}

void ImplicitSphereSource::setNoisePhaseOffset(float offset)
{
	m_noiseParamStruct->noiseFlow = offset;
}

void ImplicitSphereSource::setVelocity(vector3 velocity)
{
	m_velocity = velocity;
}

//
// private methods
//

void ImplicitSphereSource::setVeluGrid(Grid3D<float>& veluGrid)
{
#pragma omp parallel for    
	for (int i = 0; i < veluGrid.getXSize(); i++)
	{
		for (int j = 0; j < veluGrid.getYSize(); j++)
		{
			for (int k = 0; k < veluGrid.getZSize(); k++)
			{
				vector3 dist;
				float x = i * m_dx; float y = j * m_dx; float z = k * m_dx;
				dist.x = x - m_center.x;
				dist.y = y - m_center.y;
				dist.z = z - m_center.z;
				float length = sqrt(dist.x*dist.x + dist.y*dist.y + dist.z*dist.z);
				float nudensity = m_density;
				float nutemperature = m_temperature;

				float tol = m_radius - length;
				if (m_emptyInterior)
					tol = abs(tol);

				if ((!m_emptyInterior && tol > 0.0) || (m_emptyInterior && tol < m_featherAmount))
				{
					float velu = veluGrid(i, j, k);
					float vel = velu + m_dt * m_velocity.x;
					veluGrid.set(vel, i, j, k);
				}
			}
		}
	}
}

void  ImplicitSphereSource::setVelvGrid(Grid3D<float>& velvGrid)
{
#pragma omp parallel for    
	for (int i = 0; i < velvGrid.getXSize(); i++)
	{
		for (int j = 0; j < velvGrid.getYSize(); j++)
		{
			for (int k = 0; k < velvGrid.getZSize(); k++)
			{
				vector3 dist;
				float x = i * m_dx; float y = j * m_dx; float z = k * m_dx;
				dist.x = x - m_center.x;
				dist.y = y - m_center.y;
				dist.z = z - m_center.z;
				float length = sqrt(dist.x*dist.x + dist.y*dist.y + dist.z*dist.z);
				float nudensity = m_density;
				float nutemperature = m_temperature;

				float tol = m_radius - length;
				if (m_emptyInterior)
					tol = abs(tol);

				if ((!m_emptyInterior && tol > 0.0) || (m_emptyInterior && tol < m_featherAmount))
				{
					float velv = velvGrid(i, j, k);
					float vel = velv + m_dt * m_velocity.y;
					velvGrid.set(vel, i, j, k);
				}
			}
		}
	}
}

void  ImplicitSphereSource::setVelwGrid(Grid3D<float>& velwGrid)
{
#pragma omp parallel for    
	for (int i = 0; i < velwGrid.getXSize(); i++)
	{
		for (int j = 0; j < velwGrid.getYSize(); j++)
		{
			for (int k = 0; k < velwGrid.getZSize(); k++)
			{
				vector3 dist;
				float x = i * m_dx; float y = j * m_dx; float z = k * m_dx;
				dist.x = x - m_center.x;
				dist.y = y - m_center.y;
				dist.z = z - m_center.z;
				float length = sqrt(dist.x*dist.x + dist.y*dist.y + dist.z*dist.z);
				float nudensity = m_density;
				float nutemperature = m_temperature;

				float tol = m_radius - length;
				if (m_emptyInterior)
					tol = abs(tol);

				if ((!m_emptyInterior && tol > 0.0) || (m_emptyInterior && tol < m_featherAmount))
				{
					float velw = velwGrid(i, j, k);
					float vel = velw + m_dt * m_velocity.z;
					velwGrid.set(vel, i, j, k);
				}
			}
		}
	}
}