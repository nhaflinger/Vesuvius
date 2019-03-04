//
// FluidSource.h 
// Copyright (c) 2016-2018
// author: Douglas Creel
//

#ifndef FLUIDSOURCE_H
#define FLUIDSOURCE_H


#include "Vesuvius.h" 
#include "Grid3D.h"

class FluidSource
{
public:
	FluidSource();

	~FluidSource();

	virtual void update(Grid3D<float>& smokeGrid, Grid3D<float>& temperatureGrid, Grid3D<float>& fuelGrid, Grid3D<float>& veluGrid, Grid3D<float>& velvGrid, Grid3D<float>& velwGrid) = 0;

	bool isActive();

	void setUpdateBehavior(UPDATEBEHAVIOR behavior);

protected:
	void setType(FLUIDSOURCETYPE type);

	bool m_active = true;

	UPDATEBEHAVIOR m_updateBehavior = ADD;

	FLUIDSOURCETYPE m_type;

private:

};

#endif