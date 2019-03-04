//
// Vesuvius.h 
// Copyright (c) 2016-2018
// author: Douglas Creel
//

#ifndef VESUVIUS_H
#define VESUVIUS_H

#include <iostream>
#include <vector>
#include <cassert>

#include <Eigen/Sparse>

#include <omp.h>

#include <cmath>
#include <sys/timeb.h>
#include <time.h>
#include "Timer.h"
#include <openvdb/openvdb.h>
#include <openvdb/tools/LevelSetSphere.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/GridTransformer.h>
#include "Utilities.h"


using namespace Eigen;


enum MATHOPERATION {
	SETVALUE,
	ADDVALUE,
	MULTVALUE
};

typedef struct {
	MATHOPERATION setMultAdd;
	float noiseAmpScale;
	vector3 noiseFreqScale;
	float noiseFlow;
	float lacunarity;
	float gain;
	float octaves;
	float flow;
} NoiseParams;

typedef struct {
	vector3 scale;
	vector3 translate;
	vector3 rotate;
} Xform;

typedef struct {
	int xmin;
	int xmax;
	int ymin;
	int ymax;
	int zmin;
	int zmax;
} ContainerBounds;

enum LAYER {
	SOLID,
	DENSITY,
	EMPTY
};

enum ADVECTMODE {
	TRACE,
	BFECC,
	MMAC
};

enum TRACEMODE {
	SINGLESTEP,
	RUNGEKUTTA2,
	RUNGEKUTTA4
};

enum FLUIDSOURCETYPE {
	IMPLICITSPHERE,
	GEOMETRY,
	PARTICLES
};

enum UPDATEBEHAVIOR {
	COPY,
	ADD,
	SUBTRACT,
	MULTIPLY,
	DIVIDE,
	MIN,
	MAX,
	AVERAGE
};

enum BOUNDARYTYPE {
	CLOSED_BOUNDARY,
	OPEN_BOUNDARY
};

enum PRESSURESOLVER {
	PCG,
	MULTIGRID
};

enum INTERPOLATIONMODE {
	TRILINEAR,
	TRICUBIC
};

#endif
