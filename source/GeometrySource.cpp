//
// GeometrySource.cpp
// Copyright (c) 2016-2018
// author: Douglas Creel
//

#include "GeometrySource.h"

using namespace Utilities;

//
// public methods
//


GeometrySource::GeometrySource(openvdb::FloatGrid::Ptr gridPtr, bool activate, UPDATEBEHAVIOR behavior, float featherAmount, bool emptyInterior, float dx, float dt, float density, Grid3D<float>& smokeGrid, float densityScale, float temperature, Grid3D<float>& temperatureGrid, float temperatureScale, float fuel, Grid3D<float>& fuelGrid, float fuelScale, bool enableNoise, NoiseParams* noiseParam, SimplexNoise noise, Xform transform, vector3 velocity, Grid3D<float>& veluGrid, Grid3D<float>& velvGrid, Grid3D<float>& velwGrid)
{
	m_active = activate;
	m_updateBehavior = behavior;
	m_density = density;
	m_temperature = temperature;
	m_fuel = fuel;
	m_dx = dx;
	m_dt = dt;
	m_enableNoise = enableNoise;
	m_emptyInterior = emptyInterior;
	m_featherAmount = featherAmount;
	m_noiseParamStruct = noiseParam;
	m_noise = noise;
	m_gridPtr = gridPtr;
	m_scale = transform.scale;
	m_translation = transform.translate;
	m_rotation = transform.rotate;
	m_velocity = velocity;
	m_densityScale = densityScale;
	m_temperatureScale = temperatureScale;
	m_fuelScale = fuelScale;

	/*
	openvdb::math::Transform xform;
	xform = gridPtr->transform();

	// Convert the level set to a narrow-band fog volume, in which
	// interior voxels have value 1, exterior voxels have value 0, and
	// narrow-band voxels have values varying linearly from 0 to 1.

	const float outside = gridPtr->background();
	const float width = 2.0 * outside;

	// Visit and update all of the grid's active values, which correspond to
	// voxels on the narrow band.
	for (openvdb::FloatGrid::ValueOnIter iter = gridPtr->beginValueOn(); iter; ++iter)
	{
	    float dist = iter.getValue();
	    iter.setValue((outside - dist) / width);
	}

	// Visit all of the grid's inactive tile and voxel values and update the values
	// that correspond to the interior region.
	for (openvdb::FloatGrid::ValueOffIter iter = gridPtr->beginValueOff(); iter; ++iter)
	{
	    if (iter.getValue() < 0.0) 
		{
		    iter.setValue(density);
	        iter.setValueOff();
	    }
	}
	*/

#pragma omp parallel for    
	for (int i = 0; i < smokeGrid.getXSize(); i++)
	{
		for (int j = 0; j < smokeGrid.getYSize(); j++)
		{
			for (int k = 0; k < smokeGrid.getZSize(); k++)
			{
				float x = i * dx; float y = j * dx; float z = k * dx;
				float xi, yj, zk;
				xi = x * m_scale.x;
				yj = y * m_scale.y;
				zk = z * m_scale.z;
				xi = xi + m_translation.x;
				yj = yj + m_translation.y; 
				zk = zk + m_translation.z;

				const openvdb::Vec3R ijk(xi, yj, zk);

				float nudensity = m_density;
				float nutemperature = m_temperature;
				float nufuel = m_fuel;

				// Compute the value via cubic interpolation.
				openvdb::FloatGrid::ValueType v1;
				openvdb::tools::QuadraticSampler::sample(gridPtr->tree(), ijk, v1);

				float outside = m_gridPtr->background();
				float width = 2.0 * outside;
				float dist = v1;
				float v2 = (outside - dist) / width;

				if (v1 < 0.0)  
				{
					if (m_enableNoise)
					{
						float pnt[3];
						pnt[0] = m_noiseParamStruct->noiseFreqScale.x * x;
						pnt[1] = m_noiseParamStruct->noiseFreqScale.y * y;
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

GeometrySource::~GeometrySource()
{
}

void GeometrySource::update(Grid3D<float>& smokeGrid, Grid3D<float>& temperatureGrid, Grid3D<float>& fuelGrid, Grid3D<float>& veluGrid, Grid3D<float>& velvGrid, Grid3D<float>& velwGrid)
{
#pragma omp parallel for    
	for (int i = 0; i < smokeGrid.getXSize(); i++)
	{
		for (int j = 0; j < smokeGrid.getYSize(); j++)
		{
			for (int k = 0; k < smokeGrid.getZSize(); k++)
			{
				float x = i * m_dx; float y = j * m_dx; float z = k * m_dx;
				float xi, yj, zk;
				xi = x * m_scale.x;
				yj = y * m_scale.y;
				zk = z * m_scale.z;
				xi = xi + m_translation.x;
				yj = yj + m_translation.y;
				zk = zk + m_translation.z;

				const openvdb::Vec3R ijk(xi, yj, zk);

				float nudensity = m_density;
				float nutemperature = m_temperature;
				float nufuel = m_fuel;

				// Compute the value via cubic interpolation.
				openvdb::FloatGrid::ValueType v1;
				openvdb::tools::QuadraticSampler::sample(m_gridPtr->tree(), ijk, v1);

				float outside = m_gridPtr->background();
				float width = 2.0 * outside;
				float dist = v1;
				float v2 = (outside - dist) / width;

				if (v1 < 0.0)
				{
					nudensity *= v2;
					nutemperature *= v2;
					nufuel *= v2;

					float smoke = smokeGrid(i, j, k);
					float temp = temperatureGrid(i, j, k);
					float fuel = fuelGrid(i, j, k);

					if (m_enableNoise)
					{
						float pnt[3];
						pnt[0] = m_noiseParamStruct->noiseFreqScale.x * x;
						pnt[1] = m_noiseParamStruct->noiseFreqScale.y * y;
						pnt[2] = m_noiseParamStruct->noiseFreqScale.z * z;
						float noise = m_noiseParamStruct->noiseAmpScale * m_noise.turbulence(pnt, m_noiseParamStruct->lacunarity, m_noiseParamStruct->gain, m_noiseParamStruct->octaves, m_noiseParamStruct->noiseFlow);

						// add noise
						if (m_noiseParamStruct->setMultAdd == MULTVALUE)
						{
							nudensity *= noise;
							nutemperature *= noise;
							nufuel *= noise;
						}
						else if(m_noiseParamStruct->setMultAdd == ADDVALUE)
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

void GeometrySource::setDensity(float density)
{
	m_density = density;
}

void GeometrySource::setTemperature(float temperature)
{
	m_temperature = temperature;
}

void GeometrySource::setFuel(float fuel)
{
	m_fuel = fuel;
}

void GeometrySource::setDensityScale(float scale)
{
	m_densityScale = scale;
}

void GeometrySource::setTemperatureScale(float scale)
{
	m_temperatureScale = scale;
}

void GeometrySource::setFuelScale(float scale)
{
	m_fuelScale = scale;
}

void GeometrySource::setTimeStep(float step)
{
	m_dt = step;
}

void GeometrySource::setVoxelSize(float voxelsize)
{
	m_dx = voxelsize;
}

void GeometrySource::setFeatherAmount(float featherAmount)
{
	m_featherAmount = featherAmount;
}

void GeometrySource::enableNoise(bool enable)
{
	m_enableNoise = enable;
}

void GeometrySource::setNoiseAmpScale(float scale)
{
	m_noiseParamStruct->noiseAmpScale = scale;
}

void GeometrySource::setNoiseFreqScale(vector3 scale)
{
	m_noiseParamStruct->noiseFreqScale = scale;
}

void GeometrySource::setNoisePhaseOffset(float offset)
{
	m_noiseParamStruct->noiseFlow = offset;
}

void GeometrySource::setTranslation(vector3 translate)
{
	m_translation = translate;
}

void GeometrySource::setRotation(vector3 rotate)
{
	m_rotation = rotate;
}

void GeometrySource::setScale(vector3 scale)
{
	m_scale = scale;
}

void GeometrySource::setPivot(vector3 pivot)
{
	m_pivot = pivot;
}

void GeometrySource::setVelocity(vector3 velocity)
{
	m_velocity = velocity;
}

//
// private methods
//

void GeometrySource::setVeluGrid(Grid3D<float>& veluGrid)
{
#pragma omp parallel for    
	for (int i = 0; i < veluGrid.getXSize(); i++)
	{
		for (int j = 0; j < veluGrid.getYSize(); j++)
		{
			for (int k = 0; k < veluGrid.getZSize(); k++)
			{
				float x = i * m_dx; float y = j * m_dx; float z = k * m_dx;
				float xi, yj, zk;
				xi = x * m_scale.x;
				yj = y * m_scale.y;
				zk = z * m_scale.z;
				xi = xi + m_translation.x;
				yj = yj + m_translation.y;
				zk = zk + m_translation.z;

				const openvdb::Vec3R ijk(xi, yj, zk);

				float nudensity = m_density;
				float nutemperature = m_temperature;

				// Compute the value via cubic interpolation.
				openvdb::FloatGrid::ValueType v1;
				openvdb::tools::QuadraticSampler::sample(m_gridPtr->tree(), ijk, v1);

				float outside = m_gridPtr->background();
				float width = 2.0 * outside;
				float dist = v1;
				float v2 = (outside - dist) / width;

				if (v1 < 0.0)
				{
					veluGrid.set(m_velocity.x, i, j, k);
				}
			}
		}
	}
}

void  GeometrySource::setVelvGrid(Grid3D<float>& velvGrid)
{
#pragma omp parallel for    
	for (int i = 0; i < velvGrid.getXSize(); i++)
	{
		for (int j = 0; j < velvGrid.getYSize(); j++)
		{
			for (int k = 0; k < velvGrid.getZSize(); k++)
			{
				float x = i * m_dx; float y = j * m_dx; float z = k * m_dx;
				float xi, yj, zk;
				xi = x * m_scale.x;
				yj = y * m_scale.y;
				zk = z * m_scale.z;
				xi = xi + m_translation.x;
				yj = yj + m_translation.y;
				zk = zk + m_translation.z;

				const openvdb::Vec3R ijk(xi, yj, zk);

				float nudensity = m_density;
				float nutemperature = m_temperature;

				// Compute the value via cubic interpolation.
				openvdb::FloatGrid::ValueType v1;
				openvdb::tools::QuadraticSampler::sample(m_gridPtr->tree(), ijk, v1);

				float outside = m_gridPtr->background();
				float width = 2.0 * outside;
				float dist = v1;
				float v2 = (outside - dist) / width;

				if (v1 < 0.0)
				{
					velvGrid.set(m_velocity.y, i, j, k);
				}
			}
		}
	}
}

void  GeometrySource::setVelwGrid(Grid3D<float>& velwGrid)
{
#pragma omp parallel for    
	for (int i = 0; i < velwGrid.getXSize(); i++)
	{
		for (int j = 0; j < velwGrid.getYSize(); j++)
		{
			for (int k = 0; k < velwGrid.getZSize(); k++)
			{
				float x = i * m_dx; float y = j * m_dx; float z = k * m_dx;
				float xi, yj, zk;
				xi = x * m_scale.x;
				yj = y * m_scale.y;
				zk = z * m_scale.z;
				xi = xi + m_translation.x;
				yj = yj + m_translation.y;
				zk = zk + m_translation.z;

				const openvdb::Vec3R ijk(xi, yj, zk);

				float nudensity = m_density;
				float nutemperature = m_temperature;

				// Compute the value via cubic interpolation.
				openvdb::FloatGrid::ValueType v1;
				openvdb::tools::QuadraticSampler::sample(m_gridPtr->tree(), ijk, v1);

				float outside = m_gridPtr->background();
				float width = 2.0 * outside;
				float dist = v1;
				float v2 = (outside - dist) / width;

				if (v1 < 0.0)
				{
					velwGrid.set(m_velocity.z, i, j, k);
				}
			}
		}
	}
}