//
// FieldIO.h 
// Copyright (c) 2016-2018
// author: Douglas Creel
//

#ifndef FIELDIO_H
#define FIELDIO_H

#include "Vesuvius.h"
#include "Grid3D.h"


class FieldIO
{
public:
	FieldIO();

	~FieldIO();

	void initializeVDB();

	void saveVDBFields(const char* filename, double voxel_size, Grid3D<float> velu, Grid3D<float> velv, Grid3D<float> velw, Grid3D<float> smoke, Grid3D<float> temperature, Grid3D<float> fuel, Grid3D<float> heat, Grid3D<float> pressure, bool writeFireFields);

	openvdb::FloatGrid::Ptr loadVDBSource(const char* filename);

	void testVDB(const char* filename, float radius, openvdb::Vec3f center, float default_value, float voxel_size, unsigned half_width);

private:

};


#endif
