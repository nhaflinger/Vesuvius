//
// FluidSource.cpp
// Copyright (c) 2016-2018
// author: Douglas Creel
//


#include "FluidSource.h"


FluidSource::FluidSource()
{
}

FluidSource::~FluidSource()
{
}

bool FluidSource::isActive()
{
	return m_active;
}

void FluidSource::setType(FLUIDSOURCETYPE type)
{
	m_type = type;
}

void FluidSource::setUpdateBehavior(UPDATEBEHAVIOR behavior)
{
	m_updateBehavior = behavior;
}