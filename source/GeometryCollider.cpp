//
// GeometryCollider.cpp
// Copyright (c) 2016-2018
// author: Douglas Creel
//


#include "GeometryCollider.h"


GeometryCollider::GeometryCollider()
{
}

GeometryCollider::GeometryCollider(openvdb::FloatGrid::Ptr gridPtr, bool activate, float dx, Grid3D<LAYER>& layer, Grid3D<float>& volfrac, Grid3D<float>& solidu, Grid3D<float>& solidv, Grid3D<float>& solidw, Xform transform)
{
	// create geometry collision object
	m_active = activate;
	m_dx = dx;
	m_gridPtr = gridPtr;
	m_scale = transform.scale;
	m_translation = transform.translate;
	m_rotation = transform.rotate;

#pragma omp parallel for    
	for (int i = 0; i < layer.getXSize(); i++)
	{
		for (int j = 0; j < layer.getYSize(); j++)
		{
			for (int k = 0; k < layer.getZSize(); k++)
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

				// Compute the value via cubic interpolation.
				openvdb::FloatGrid::ValueType v1;
				openvdb::tools::QuadraticSampler::sample(gridPtr->tree(), ijk, v1);

				float outside = m_gridPtr->background();
				float width = 2.0 * outside;
				float dist = v1;
				float v2 = (outside - dist) / width;

				if (v1 < 0.0)
				{
					layer.set(SOLID, i, j, k);
				}
			}
		}
	}

	// compute volume fraction of air in grid cell
#pragma omp parallel for    
	for (int i = 0; i < volfrac.getXSize(); i++)
	{
		for (int j = 0; j < volfrac.getYSize(); j++)
		{
			for (int k = 0; k < volfrac.getZSize(); k++)
			{
				// sample using a 2x2 grid to determine volume fraction of air in grid cell

				float x = i * dx; float y = j * dx; float z = k * dx;
				float xi, yj, zk;
				xi = x * m_scale.x;
				yj = y * m_scale.y;
				zk = z * m_scale.z;
				xi = xi + m_translation.x;
				yj = yj + m_translation.y;
				zk = zk + m_translation.z;

				const openvdb::Vec3R ijk(xi, yj, zk);

				// Compute the value via cubic interpolation.
				openvdb::FloatGrid::ValueType v1;
				openvdb::tools::QuadraticSampler::sample(gridPtr->tree(), ijk, v1);

				float outside = m_gridPtr->background();
				float width = 2.0 * outside;
				float dist = v1;
				float v2 = (outside - dist) / width;

				if (v1 < 0.0)
				{
					volfrac.set(SOLID, i, j, k);
				}
			}
		}
	}
}

GeometryCollider::~GeometryCollider()
{
}
