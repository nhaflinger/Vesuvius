//
// ImplicitSphereCollider.cpp
// Copyright (c) 2016-2018
// author: Douglas Creel
//


#include "ImplicitSphereCollider.h"


ImplicitSphereCollider::ImplicitSphereCollider()
{
}

ImplicitSphereCollider::ImplicitSphereCollider(bool activate, float dx, vector3 center, float radius, Grid3D<LAYER>& layer, Grid3D<float>& volfrac, Grid3D<float>& solidu, Grid3D<float>& solidv, Grid3D<float>& solidw)
{
	// create implicit sphere collision object
	m_active = activate;
	m_center.x = center.x; m_center.y = center.y; m_center.z = center.z;
	m_radius = radius;
	m_dx = dx;

#pragma omp parallel for    
	for (int i = 0; i < layer.getXSize(); i++)
	{
		for (int j = 0; j < layer.getYSize(); j++)
		{
			for (int k = 0; k < layer.getZSize(); k++)
			{
				vector3 dist;
				float x = i * dx; float y = j * dx; float z = k * dx;
				dist.x = x - center.x;
				dist.y = y - center.y;
				dist.z = z - center.z;
				float length = sqrt(dist.x*dist.x + dist.y*dist.y + dist.z*dist.z);
				float tol = radius - length;
				if (tol > 0.5f * dx)
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
			}
		}
	}
}

ImplicitSphereCollider::~ImplicitSphereCollider()
{
}