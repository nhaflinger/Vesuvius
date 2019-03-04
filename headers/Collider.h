//
// Collider.h 
// Copyright (c) 2016-2018
// author: Douglas Creel
//

#ifndef COLLIDER_H
#define COLLIDER_H


#include "Vesuvius.h" 
#include "Grid3D.h"

class Collider
{
public:
	Collider();

	~Collider();

	bool isActive();

protected:
	bool m_active = true;

private:

};

#endif
