//
// GeometryCollider.h 
// Copyright (c) 2016-2018
// author: Douglas Creel
//

#ifndef GEOMETRYCOLLIDER_H
#define GEOMETRYCOLLIDER_H


#include "Vesuvius.h" 
#include "Grid3D.h"
#include "Collider.h"


class GeometryCollider : public Collider
{
public:
	GeometryCollider();

	GeometryCollider(openvdb::FloatGrid::Ptr gridPtr, bool activate, float dx, Grid3D<LAYER>& layer, Grid3D<float>& volfrac, Grid3D<float>& solidu, Grid3D<float>& solidv, Grid3D<float>& solidw, Xform transform);

	~GeometryCollider();


	inline void setTimeStep(float step)
	{
		m_dt = step;
	}

	inline void setVoxelSize(float voxelsize)
	{
		m_dx = voxelsize;
	}

	inline void setTranslation(vector3 translate)
	{
		m_translation = translate;
	}

	inline void setRotation(vector3 rotate)
	{
		m_rotation = rotate;
	}

	inline void setScale(vector3 scale)
	{
		m_scale = scale;
	}

	inline void setPivot(vector3 pivot)
	{
		m_pivot = pivot;
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
	double m_dx;

	double m_dt;

	openvdb::FloatGrid::Ptr m_gridPtr;

	vector3 m_translation = { 0, 0, 0 };

	vector3 m_rotation = { 0, 0, 0 };

	vector3 m_scale = { 1, 1, 1 };

	vector3 m_pivot = { 0, 0, 0 };

	bool m_staticCollider = true;
};

#endif