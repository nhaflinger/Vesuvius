//
// FieldIO.cpp
// Copyright (c) 2016-2018
// author: Douglas Creel
//

#include "FieldIO.h"


FieldIO::FieldIO()
{
}

FieldIO::~FieldIO()
{
}

void FieldIO::initializeVDB()
{
	openvdb::initialize();
}

openvdb::FloatGrid::Ptr FieldIO::loadVDBSource(const char* filename)
{
	openvdb::io::File file(filename);

	file.open();

	openvdb::GridBase::Ptr baseGrid;
	for (openvdb::io::File::NameIterator nameIter = file.beginName(); nameIter != file.endName(); ++nameIter)
	{
		if (nameIter.gridName() == "surface") 
		{
			baseGrid = file.readGrid(nameIter.gridName());
		}
		//else
		//{
			//std::cout << "other grid names: " << nameIter.gridName() << std::endl;
		//}
	}

	file.close();
	
	openvdb::FloatGrid::Ptr grid = openvdb::gridPtrCast<openvdb::FloatGrid>(baseGrid);

	return grid;
}

void FieldIO::saveVDBFields(const char* filename, double voxel_size, Grid3D<float> velu, Grid3D<float> velv, Grid3D<float> velw, Grid3D<float> smoke, Grid3D<float> temperature, Grid3D<float> fuel, Grid3D<float> heat, Grid3D<float> pressure, bool writeFireFields)
{
	double threshold = 1.0e-06;

	// populate smokeGrid with values	
	openvdb::FloatGrid::Ptr smokeGrid = openvdb::FloatGrid::create(0.0);
	smokeGrid->setName("smoke");
	smokeGrid->setTransform(openvdb::math::Transform::createLinearTransform(voxel_size));
	smokeGrid->setGridClass(openvdb::GRID_FOG_VOLUME);

	openvdb::FloatGrid::Accessor accessor2 = smokeGrid->getAccessor();
	for (int i = 0; i < smoke.getXSize(); i++)
	{
		for (int j = 0; j < smoke.getYSize(); j++)
		{
			for (int k = 0; k < smoke.getZSize(); k++)
			{
				openvdb::Coord xyz(i, j, k);
				if (smoke(i, j, k) > threshold)
					accessor2.setValue(xyz, smoke(i, j, k));
			}
		}
	}

	// populate temperatureGrid with values
	openvdb::FloatGrid::Ptr temperatureGrid = openvdb::FloatGrid::create(0.0);
	temperatureGrid->setName("temperature");
	temperatureGrid->setTransform(openvdb::math::Transform::createLinearTransform(voxel_size));
	temperatureGrid->setGridClass(openvdb::GRID_FOG_VOLUME);

	openvdb::FloatGrid::Accessor accessor3 = temperatureGrid->getAccessor();
	for (int i = 0; i < temperature.getXSize(); i++)
	{
		for (int j = 0; j < temperature.getYSize(); j++)
		{
			for (int k = 0; k < temperature.getZSize(); k++)
			{
				openvdb::Coord xyz(i, j, k);
				if (temperature(i, j, k) > threshold)
					accessor3.setValue(xyz, temperature(i, j, k));
			}
		}
	}

	// populate fuelGrid with values	
	openvdb::FloatGrid::Ptr fuelGrid = openvdb::FloatGrid::create(0.0);
	fuelGrid->setName("fuel");
	fuelGrid->setTransform(openvdb::math::Transform::createLinearTransform(voxel_size));
	fuelGrid->setGridClass(openvdb::GRID_FOG_VOLUME);

	openvdb::FloatGrid::Accessor accessor6 = fuelGrid->getAccessor();
	for (int i = 0; i < fuel.getXSize(); i++)
	{
		for (int j = 0; j < fuel.getYSize(); j++)
		{
			for (int k = 0; k < fuel.getZSize(); k++)
			{
				openvdb::Coord xyz(i, j, k);
				if (fuel(i, j, k) > threshold)
					accessor6.setValue(xyz, fuel(i, j, k));
			}
		}
	}

	// populate heatGrid with values	
	openvdb::FloatGrid::Ptr heatGrid = openvdb::FloatGrid::create(0.0);
	heatGrid->setName("heat");
	heatGrid->setTransform(openvdb::math::Transform::createLinearTransform(voxel_size));
	heatGrid->setGridClass(openvdb::GRID_FOG_VOLUME);

	openvdb::FloatGrid::Accessor accessor8 = heatGrid->getAccessor();
	for (int i = 0; i < heat.getXSize(); i++)
	{
		for (int j = 0; j < heat.getYSize(); j++)
		{
			for (int k = 0; k < heat.getZSize(); k++)
			{
				openvdb::Coord xyz(i, j, k);
				if (heat(i, j, k) > threshold)
					accessor8.setValue(xyz, heat(i, j, k));
			}
		}
	}
/*
	// populate pressureGrid with values	
	openvdb::FloatGrid::Ptr pressureGrid = openvdb::FloatGrid::create(0.0);
	pressureGrid->setName("pressure");
	pressureGrid->setTransform(openvdb::math::Transform::createLinearTransform(voxel_size));
	pressureGrid->setGridClass(openvdb::GRID_FOG_VOLUME);

	openvdb::FloatGrid::Accessor accessor7 = pressureGrid->getAccessor();
	for (int i = 0; i < pressure.getXSize(); i++)
	{
		for (int j = 0; j < pressure.getYSize(); j++)
		{
			for (int k = 0; k < pressure.getZSize(); k++)
			{
				openvdb::Coord xyz(i, j, k);
				if (abs(pressure(i, j, k)) > threshold)
					accessor7.setValue(xyz, pressure(i, j, k));
			}
		}
	}
*/
	// populate velocity grids with values
	openvdb::FloatGrid::Ptr velxGrid = openvdb::FloatGrid::create(0.0);
	velxGrid->setName("vel.x");
	velxGrid->setTransform(openvdb::math::Transform::createLinearTransform(voxel_size));
	velxGrid->setGridClass(openvdb::GRID_STAGGERED);

	openvdb::FloatGrid::Accessor accessor1 = velxGrid->getAccessor();
	for (int i = 0; i < velu.getXSize(); i++)
	{
		for (int j = 0; j < velu.getYSize(); j++)
		{
			for (int k = 0; k < velu.getZSize(); k++)
			{
				openvdb::Coord xyz(i, j, k);
				if ((velu(i, j, k)) > threshold)
					accessor1.setValue(xyz, velu(i, j, k));
			}
		}
	}

	openvdb::FloatGrid::Ptr velyGrid = openvdb::FloatGrid::create(0.0);
	velyGrid->setName("vel.y");
	velyGrid->setTransform(openvdb::math::Transform::createLinearTransform(voxel_size));
	velyGrid->setGridClass(openvdb::GRID_STAGGERED);

	openvdb::FloatGrid::Accessor accessor4 = velyGrid->getAccessor();
	for (int i = 0; i < velv.getXSize(); i++)
	{
		for (int j = 0; j < velv.getYSize(); j++)
		{
			for (int k = 0; k < velv.getZSize(); k++)
			{
				openvdb::Coord xyz(i, j, k);
				if (abs(velv(i, j, k)) > threshold)
					accessor4.setValue(xyz, velv(i, j, k));
			}
		}
	}

	openvdb::FloatGrid::Ptr velzGrid = openvdb::FloatGrid::create(0.0);
	velzGrid->setName("vel.z");
	velzGrid->setTransform(openvdb::math::Transform::createLinearTransform(voxel_size));
	velzGrid->setGridClass(openvdb::GRID_STAGGERED);

	openvdb::FloatGrid::Accessor accessor5 = velzGrid->getAccessor();
	for (int i = 0; i < velw.getXSize(); i++)
	{
		for (int j = 0; j < velw.getYSize(); j++)
		{
			for (int k = 0; k < velw.getZSize(); k++)
			{
				openvdb::Coord xyz(i, j, k);
				if (abs(velw(i, j, k)) > threshold)
					accessor5.setValue(xyz, velw(i, j, k));
			}
		}
	}

	// Add the grid pointers to a container.
	openvdb::GridPtrVec grids;
	grids.push_back(smokeGrid);
	grids.push_back(temperatureGrid);
	grids.push_back(fuelGrid);
	grids.push_back(heatGrid);
	//grids.push_back(pressureGrid);
	grids.push_back(velxGrid);
	grids.push_back(velyGrid);
	grids.push_back(velzGrid);

	// Write out the contents of the container
	openvdb::io::File file(filename);
	file.write(grids);
	file.close();
}

void FieldIO::testVDB(const char* filename, float radius, openvdb::Vec3f center, float default_value, float voxel_size, unsigned half_width)
{
	// Create a FloatGrid and populate it with a narrow-band
	// signed distance field of a sphere.
	openvdb::FloatGrid::Ptr grid =
		openvdb::tools::createLevelSetSphere<openvdb::FloatGrid>(
			/*radius=*/radius, /*center=*/center,
			/*voxel size=*/voxel_size, /*width=*/half_width);

	// Associate some metadata with the grid.
	grid->insertMeta("radius", openvdb::FloatMetadata(radius));

	// Name the grid "LevelSetSphere".
	grid->setName("LevelSetSphere");

	// Create a VDB file object.
	openvdb::io::File file(filename);

	// Add the grid pointer to a container.
	openvdb::GridPtrVec grids;
	grids.push_back(grid);

	// Write out the contents of the container.
	file.write(grids);
	file.close();
}

