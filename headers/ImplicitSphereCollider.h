//
// ImplicitSphereCollider.h 
// Copyright (c) 2016-2018
// author: Douglas Creel
//

#ifndef IMPLICITSPHERECOLLIDER_H
#define IMPLICITSPHERECOLLIDER_H


#include "Vesuvius.h" 
#include "Grid3D.h"
#include "Collider.h"


class ImplicitSphereCollider : public Collider
{
public:
	ImplicitSphereCollider();

	ImplicitSphereCollider(bool activate, float dx, vector3 center, float radius, Grid3D<LAYER>& layer, Grid3D<float>& volfrac, Grid3D<float>& solidu, Grid3D<float>& solidv, Grid3D<float>& solidw);

	~ImplicitSphereCollider();

	inline void setCenter(vector3 center)
	{
		m_center = center;
	}

	inline void setRadius(double radius)
	{
		m_radius = radius;
	}

	inline void setTimeStep(double step)
	{
		m_dt = step;
	}

	inline void setVoxelSize(double voxelsize)
	{
		m_dx = voxelsize;
	}

	inline bool setStaticFlag(bool flag)
	{
		m_staticCollider = flag;
	}

	inline bool getStaticFlag()
	{
		return m_staticCollider;
	}

private:
	float m_dx;

	float m_dt;

	vector3 m_center;

	float m_radius;

	bool m_staticCollider = true;

};
#endif