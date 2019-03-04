//
// GasSolver.cpp
// Copyright (c) 2016-2018
// author: Douglas Creel
//

#include "GasSolver.h"
#include "FieldIO.h"

#define FIELD_THRESHOLD 1.0e-06
#define VELOCITY_THRESHOLD 1.0e-16


GasSolver::GasSolver()
{
}

GasSolver::GasSolver(int xsize, int ysize, int zsize, float dx)
{
	m_xsize = xsize;
	m_ysize = ysize;
	m_zsize = zsize;
	m_dx = dx;
	m_simBounds = { 0, m_xsize - 1, 0, m_ysize - 1, 0, m_zsize - 1 };
}

GasSolver::~GasSolver()
{
}

//
// Public methods
//

void GasSolver::initialize()
{
	if (!m_initialized)
	{
		m_initializeSimulation();
	}
}

void GasSolver::addBodyForce(vector3 force)
{
	m_gravity = force;
}

void GasSolver::update()
{
	int maxsteps = clamp(m_maxSubSteps, 1, 100);

	if (m_simTime == 0)
	{
		setSolveRegion();

		// initialize pressure coefficients
		computePressureCoefficients();
	}

	for (int i = 0; i < m_subSteps; i++)
	{
		float dt = computeCFL();
		setTimeStep(getSubSteps());
		float timestep = m_timeStep;
		int substeps = 1;
		while (timestep >= dt && substeps < maxsteps)
		{
			timestep -= dt;
			substeps += 1;
		}
		m_timeStep /= substeps;
		for (int j = 0; j < substeps; j++)
		{
			advanceOneStep(); 
			m_simTime += m_timeStep;
		}
	}

	setTimeStep(getSubSteps());
}

void GasSolver::advanceOneStep()
{
	// update density advection region
	setSolveRegion();

	// resize container?
	if (m_autoResizeContainer)
		autoResizeContainer();

	// update collision objects
#pragma omp parallel sections
	{
#pragma omp section
		updateColliders();

		// update fluid sources
#pragma omp section
		updateFluidSources();
	}

	// apply forces
	applyBouyancyForce();

	// apply wind force
	if (m_enableWindForce)
		applyWindForce();

	// apply turbulence force
	if (m_enableTurbulence)
		applyTurbulence();

	// apply vorticity confinement force
	if (m_enableVorticityConfinement && m_vorticityScale != 0)
		applyVorticityConfinement();

	// apply shredding to fire
	if (m_enableShredding && m_shreddingAmount > 0)
		applyShredding();

	// apply viscosity
	if (m_enableViscosity && m_kinematicViscosity != 0)
		applyViscosity();

	// advect velocity 
	advectVelocity();

	// extrapolate velocities
	applyBoundaryConditions();

	// pressure projection
	project();

	// extrapolate velocities
	applyBoundaryConditions();

	// advect scalar quantities (smoke, temperature, fire)
	advectScalars();

	// apply diffusion
	if (m_diffusion != 0)
		applyDiffusion();

	// apply dissipation
	if (m_enableDissipation && m_dissipation != 0)
		applyDissipation();

	// apply cooling
	if (m_cooling > 0)
		applyCooling();

	// apply fuel dissipation
	if (m_enableFuelDissipation && m_fuelDissipation != 0)
		applyFuelDissipation();

	// apply fuel cooling
	if (m_flameCooling != 0)
		applyFlameCooling();
}

void GasSolver::saveState(const char* filename)
{
	// open VDB file for writing
	FieldIO inputOutput;
	inputOutput.initializeVDB();
	inputOutput.saveVDBFields(filename, m_dx, m_velu, m_velv, m_velw, m_smoke, m_temperature, m_fuel, m_heat, m_pressure, m_enableFire);
}

void GasSolver::implicitSphereSource(UPDATEBEHAVIOR behavior, bool emptyInterior, float featherAmount, vector3 center, float radius, float density, float densityScale, float temperature, float temperatureScale, float fuel, float fuelScale, bool enableNoise, NoiseParams* noiseParam, SimplexNoise noise, vector3 velocity)
{
	// clamp scale factors
	if (!densityScale) densityScale = FIELD_THRESHOLD;
	if (!temperatureScale) temperatureScale = FIELD_THRESHOLD;
	if (!fuelScale) fuelScale = FIELD_THRESHOLD;

	// create implicit sphere with a set smoke density and temperature
	ImplicitSphereSource smokeSource(true, behavior, featherAmount, emptyInterior, m_dx, m_timeStep, center, radius, density, m_smoke, densityScale, temperature, m_temperature, temperatureScale, fuel, m_fuel, fuelScale, enableNoise, noiseParam, noise, velocity, m_velu,  m_velv, m_velw);

	m_implicitSphereSources.push_back(smokeSource);
}

void GasSolver::implicitSphereCollider(vector3 center, float radius)
{
	// create implicit sphere collider
	ImplicitSphereCollider collider(true, m_dx, center, radius, m_layer, m_volfrac, m_solidu, m_solidv, m_solidw);

	m_implicitSphereColliders.push_back(collider);
}

void GasSolver::geometrySource(UPDATEBEHAVIOR behavior, bool emptyInterior, float featherAmount, const char* geometryFile, float density, float densityScale, float temperature, float temperatureScale, float fuel, float fuelScale, bool enableNoise, NoiseParams* noiseparam, SimplexNoise noise, Xform transform, vector3 velocity)
{
	// open VDB file for reading
	FieldIO inputOutput;
	inputOutput.initializeVDB();

	openvdb::FloatGrid::Ptr gridPtr = inputOutput.loadVDBSource(geometryFile);

	// clamp scale factors
	if (!densityScale) densityScale = FIELD_THRESHOLD;
	if (!temperatureScale) temperatureScale = FIELD_THRESHOLD;
	if (!fuelScale) fuelScale = FIELD_THRESHOLD;

	GeometrySource smokeSource(gridPtr, true, behavior, emptyInterior, featherAmount, m_dx, m_timeStep, density, m_smoke, densityScale, temperature, m_temperature, temperatureScale, fuel, m_fuel, fuelScale, enableNoise, noiseparam, noise, transform, velocity, m_velu, m_velv, m_velw);

	m_geometrySources.push_back(smokeSource);
}

void GasSolver::geometryCollider(const char* geometryFile, Xform transform)
{
	// open VDB file for reading
	FieldIO inputOutput;
	inputOutput.initializeVDB();

	openvdb::FloatGrid::Ptr gridPtr = inputOutput.loadVDBSource(geometryFile);

	GeometryCollider collider(gridPtr, true, m_dx, m_layer, m_volfrac, m_solidu, m_solidv, m_solidw, transform);

	m_geometryColliders.push_back(collider);
}

float GasSolver::computeCFL()
{
	float gravity = sqrt(m_gravity.x * m_gravity.x + m_gravity.y * m_gravity.y + m_gravity.z * m_gravity.z);
	float a0 = m_dx * gravity;
	float a1 = (m_velu.maxval() * m_velu.maxval()) + (m_velv.maxval() * m_velv.maxval()) + (m_velw.maxval() * m_velw.maxval());
	float maxv2 = max(a0, a1);

	if (maxv2 < VELOCITY_THRESHOLD) maxv2 = VELOCITY_THRESHOLD;

	return m_cflScale * m_dx / sqrt(maxv2);
}

//
// Private methods
//

void GasSolver::m_initializeSimulation()
{
	int xdim = getXSize();
	int ydim = getYSize();
	int zdim = getZSize();

#pragma omp parallel sections
	{
#pragma omp section
		m_velu.init(xdim + 1, ydim, zdim, 0.f);
#pragma omp section
		m_velv.init(xdim, ydim + 1, zdim, 0.f);
#pragma omp section
		m_velw.init(xdim, ydim, zdim + 1, 0.f);

#pragma omp section
		m_smoke.init(xdim, ydim, zdim, 0.f);
#pragma omp section
		m_temperature.init(xdim, ydim, zdim, 0.f);
#pragma omp section
		m_fuel.init(xdim, ydim, zdim, 0.f);
#pragma omp section
		m_heat.init(xdim, ydim, zdim, 0.f);
#pragma omp section
		m_density.init(xdim, ydim, zdim, 0.f);
#pragma omp section
		m_layer.init(xdim, ydim, zdim, DENSITY);
#pragma omp section
		m_volfrac.init(xdim + 1, ydim + 1, zdim + 1, DENSITY);

#pragma omp section
		m_solidu.init(xdim + 1, ydim, zdim, 0.f);
#pragma omp section
		m_solidv.init(xdim, ydim + 1, zdim, 0.f);
#pragma omp section
		m_solidw.init(xdim, ydim, zdim + 1, 0.f);

#pragma omp section
		m_restu.init(xdim, ydim, zdim, 0.f);
#pragma omp section
		m_restv.init(xdim, ydim, zdim, 0.f);
#pragma omp section
		m_restw.init(xdim, ydim, zdim, 0.f);
	}

#pragma omp parallel for   
	for (int i = 0; i < xdim; i++)
	{
		for (int j = 0; j < ydim; j++)
		{
			for (int k = 0; k < zdim; k++)
			{
				vector3 pos;
				pos.x = i * m_dx; pos.y = j * m_dx; pos.z = k * m_dx;
				m_restu.set(pos.x, i, j, k);
				m_restv.set(pos.y, i, j, k);
				m_restw.set(pos.z, i, j, k);
			}
        }
    }

	m_pressure.init(xdim, ydim, zdim, 0.f);
	int nvoxels = m_pressure.getNumVoxels();

#pragma omp parallel sections
	{
#pragma omp section
		m_divergence.resize(nvoxels);
#pragma omp section
		m_coefficients.resize(nvoxels, nvoxels);
#pragma omp section
		m_guess.resize(nvoxels);
#pragma omp section
		m_guess.setZero();
	}
}

bool GasSolver::isInitialized()
{
	return m_initialized;
}

vector3 GasSolver::traceParticle(float x, float y, float z, float dt, float offset)
{
	vector3 vel = getVelocity(x, y, z, offset);
	vel = getVelocity(x + 0.5f * dt * vel.x, y + 0.5f * dt * vel.y, z + 0.5f * dt * vel.z, offset);

	vector3 pnt;
	pnt.x = x + dt * vel.x;
	pnt.y = y + dt * vel.y;
	pnt.z = z + dt * vel.z;

	// clamp to bounds of grid
	pnt.x = clamp(pnt.x, 0.0, m_velu.getXSize()*m_dx);
	pnt.y = clamp(pnt.y, 0.0, m_velv.getYSize()*m_dx);
	pnt.z = clamp(pnt.z, 0.0, m_velw.getZSize()*m_dx);

	return pnt;
}

// Back and forth error correction and compensation 
vector3 GasSolver::bfecc(float x, float y, float z, float dt, float offset)
{
	vector3 pnt, pntnew, pntback, vel, velback, velnew, velorig;

	if (m_traceMode == RUNGEKUTTA2)    // This uses Runge-Kutta 2 
	{	
		velorig = getVelocity(x, y, z, offset);
		pntnew = traceParticle(x, y, z, dt, offset);
		pntback = traceParticle(pntnew.x, pntnew.y, pntnew.z, -dt, offset);
		velback = getVelocity(pntback.x, pntback.y, pntback.z, offset);

		velnew.x = velorig.x + 0.5f * (velorig.x - velback.x);
		velnew.y = velorig.y + 0.5f * (velorig.y - velback.y);
		velnew.z = velorig.z + 0.5f * (velorig.z - velback.z);

		vel = getVelocity(x + 0.5f * dt * velnew.x, y + 0.5f * dt * velnew.y, z + 0.5f * dt * velnew.z, offset);
		pnt.x = x + dt * vel.x; pnt.y = y + dt * vel.y; pnt.z = z + dt * vel.z;
	}
	else    // This uses Euler step
	{
		vector3 velorig = getVelocity(x, y, z, offset);

		pnt.x = x + dt * velorig.x; pnt.y = y + dt * velorig.y; pnt.z = z + dt * velorig.z;
		vector3 velnew = getVelocity(pnt.x, pnt.y, pnt.z, offset);

		pnt.x = x - dt * velnew.x; pnt.y = y - dt * velnew.y; pnt.z = z - dt * velnew.z;
		vector3 velback = getVelocity(pnt.x, pnt.y, pnt.z, offset);

		velnew.x = velorig.x + 0.5f * (velorig.x - velback.x);
		velnew.y = velorig.y + 0.5f * (velorig.y - velback.y);
		velnew.z = velorig.z + 0.5f * (velorig.z - velback.z);

		pnt.x = x + dt * velnew.x; pnt.y = y + dt * velnew.y; pnt.z = z + dt * velnew.z;
	}

	// clamp to bounds of grid
	pnt.x = clamp(pnt.x, 0.0, m_velu.getXSize()*m_dx);
	pnt.y = clamp(pnt.y, 0.0, m_velv.getYSize()*m_dx);
	pnt.z = clamp(pnt.z, 0.0, m_velw.getZSize()*m_dx);

	return pnt;
}

// Modified MacCormack TBD
vector3 GasSolver::modifiedMacCormack(float x, float y, float z, float dt, float offset)
{
	vector3 vel = getVelocity(x, y, z, offset);
	vel = getVelocity(x - 0.5f * dt * vel.x, y - 0.5f * dt * vel.y, z - 0.5f * dt * vel.z, offset);

	vector3 pnt;
	pnt.x = x - dt * vel.x;
	pnt.y = y - dt * vel.y;
	pnt.z = z - dt * vel.z;

	// clamp to bounds of grid
	pnt.x = clamp(pnt.x, 0.0, m_velu.getXSize()*m_dx);
	pnt.y = clamp(pnt.y, 0.0, m_velv.getYSize()*m_dx);
	pnt.z = clamp(pnt.z, 0.0, m_velw.getZSize()*m_dx);

	return pnt;
}

vector3 GasSolver::getVelocity(float x, float y, float z, float offset)
{
	vector3 vel; 

	if (m_velocityInterpolationMode == TRICUBIC)
	{
		vel.x = m_velu.tricubicInterpolation(x / m_dx + offset, y / m_dx, z / m_dx);
		vel.y = m_velv.tricubicInterpolation(x / m_dx, y / m_dx + offset, z / m_dx);
		vel.z = m_velw.tricubicInterpolation(x / m_dx, y / m_dx, z / m_dx + offset);
	}
	else
	{
		vel.x = m_velu.trilinearInterpolation(x / m_dx + offset, y / m_dx, z / m_dx);
		vel.y = m_velv.trilinearInterpolation(x / m_dx, y / m_dx + offset, z / m_dx);
		vel.z = m_velw.trilinearInterpolation(x / m_dx, y / m_dx, z / m_dx + offset);
	}

	return vel;
}

void GasSolver::advectVelocity()
{
	Grid3D<float> velu;
	Grid3D<float> velv;
	Grid3D<float> velw;

#pragma omp parallel sections
	{
#pragma omp section
		velu.init(m_velu.getXSize(), m_velu.getYSize(), m_velu.getZSize(), 0.f);
#pragma omp section
		velv.init(m_velv.getXSize(), m_velv.getYSize(), m_velv.getZSize(), 0.f);
#pragma omp section
		velw.init(m_velw.getXSize(), m_velw.getYSize(), m_velw.getZSize(), 0.f);
	}
 
#pragma omp parallel for   
	for (int i = 0; i < m_velu.getXSize(); i++)
	{
		for (int j = 0; j < m_velu.getYSize(); j++)
		{
			for (int k = 0; k < m_velu.getZSize(); k++)
			{
				if (m_layer(i, j, k) != SOLID)
				{
					float x, y, z;
					vector3 newpnt;
					vector3 newvel;

					x = i * m_dx - 0.5f * m_dx;
					y = j * m_dx;
					z = k * m_dx;

					switch (m_advectMode)
					{
					case BFECC:
						newpnt = bfecc(x, y, z, -m_timeStep, 0.5f);
						break;
					default:
						newpnt = traceParticle(x, y, z, -m_timeStep, 0.5f);
					}

					newvel = getVelocity(newpnt.x, newpnt.y, newpnt.z, 0.5f);

					velu.set(newvel.x, i, j, k);
				}
			}
		}
	}

#pragma omp parallel for   
	for (int i = 0; i < m_velv.getXSize(); i++)
	{
		for (int j = 0; j < m_velv.getYSize(); j++)
		{
			for (int k = 0; k < m_velv.getZSize(); k++)
			{
				if (m_layer(i, j, k) != SOLID)
				{
					float x, y, z;
					vector3 newpnt;
					vector3 newvel;

					x = i * m_dx;
					y = j * m_dx - 0.5f * m_dx;
					z = k * m_dx;

					switch (m_advectMode)
					{
					case BFECC:
						newpnt = bfecc(x, y, z, -m_timeStep, 0.5f);
						break;
					default:
						newpnt = traceParticle(x, y, z, -m_timeStep, 0.5f);
					}

					newvel = getVelocity(newpnt.x, newpnt.y, newpnt.z, 0.5f);

					velv.set(newvel.y, i, j, k);
				}
			}
		}
	}

#pragma omp parallel for   
	for (int i = 0; i < m_velw.getXSize(); i++)
	{
		for (int j = 0; j < m_velw.getYSize(); j++)
		{
			for (int k = 0; k < m_velw.getZSize(); k++)
			{
				if (m_layer(i, j, k) != SOLID)
				{
					float x, y, z;
					vector3 newpnt;
					vector3 newvel;

					x = i * m_dx;
					y = j * m_dx;
					z = k * m_dx - 0.5f * m_dx;

					switch (m_advectMode)
					{
					case BFECC:
						newpnt = bfecc(x, y, z, -m_timeStep, 0.5f);
						break;
					default:
						newpnt = traceParticle(x, y, z, -m_timeStep, 0.5f);
					}

					newvel = getVelocity(newpnt.x, newpnt.y, newpnt.z, 0.5f);

					velw.set(newvel.z, i, j, k);
				}
			}
		}
	}

#pragma omp parallel for 
	for (int i = 0; i < m_velu.getNumVoxels(); i++)
	{
		m_velu.setIndexValue(i, velu[i]);
	}

#pragma omp parallel for    
	for (int i = 0; i < m_velv.getNumVoxels(); i++)
	{
		m_velv.setIndexValue(i, velv[i]);
	}

#pragma omp parallel for    
	for (int i = 0; i < m_velw.getNumVoxels(); i++)
	{
		m_velw.setIndexValue(i, velw[i]);
	}
}

float GasSolver::getSmoke(float x, float y, float z)
{
	float smoke;
	if (m_densityInterpolationMode == TRICUBIC)
	{
		smoke = m_smoke.tricubicInterpolation(x / m_dx, y / m_dx, z / m_dx);
	}
	else
	{
		smoke = m_smoke.trilinearInterpolation(x / m_dx, y / m_dx, z / m_dx);
	}
	return smoke;
}

float GasSolver::getTemperature(float x, float y, float z)
{
	float temperature;
	if (m_densityInterpolationMode == TRICUBIC)
	{
		temperature = m_temperature.tricubicInterpolation(x / m_dx, y / m_dx, z / m_dx);
	}
	else
	{
		temperature = m_temperature.trilinearInterpolation(x / m_dx, y / m_dx, z / m_dx);
	}
	return temperature;
}

float GasSolver::getFuel(float x, float y, float z)
{
	float fuel;
	if (m_densityInterpolationMode == TRICUBIC)
	{
		fuel = m_fuel.tricubicInterpolation(x / m_dx, y / m_dx, z / m_dx);
	}
	else
	{
		fuel = m_fuel.trilinearInterpolation(x / m_dx, y / m_dx, z / m_dx);
	}
	return fuel;
}

float GasSolver::getHeat(float x, float y, float z)
{
	float heat;
	if (m_densityInterpolationMode == TRICUBIC)
	{
		heat = m_heat.tricubicInterpolation(x / m_dx, y / m_dx, z / m_dx);
	}
	else
	{
		heat = m_heat.trilinearInterpolation(x / m_dx, y / m_dx, z / m_dx);
	}
	return heat;
}

vector3 GasSolver::getRest(float x, float y, float z)
{
	vector3 rest;
	if (m_densityInterpolationMode == TRICUBIC)
	{
		rest.x = m_restu.tricubicInterpolation(x / m_dx, y / m_dx, z / m_dx);
		rest.y = m_restv.tricubicInterpolation(x / m_dx, y / m_dx, z / m_dx);
		rest.z = m_restw.tricubicInterpolation(x / m_dx, y / m_dx, z / m_dx);
	}
	else
	{
		rest.x = m_restu.trilinearInterpolation(x / m_dx, y / m_dx, z / m_dx);
		rest.y = m_restv.trilinearInterpolation(x / m_dx, y / m_dx, z / m_dx);
		rest.z = m_restw.trilinearInterpolation(x / m_dx, y / m_dx, z / m_dx);
	}
	return rest;
}

void GasSolver::advectScalars()
{
	Grid3D<float> smoke;
	Grid3D<float> temperature;
	Grid3D<float> restu;
	Grid3D<float> restv;
	Grid3D<float> restw;
	Grid3D<float> fuel;
	Grid3D<float> heat;

#pragma omp parallel sections
	{
#pragma omp section
		smoke.init(m_smoke.getXSize(), m_smoke.getYSize(), m_smoke.getZSize(), 0.f);
#pragma omp section
		temperature.init(m_temperature.getXSize(), m_temperature.getYSize(), m_temperature.getZSize(), 0.f);
#pragma omp section
		restu.init(m_restu.getXSize(), m_restu.getYSize(), m_restu.getZSize(), 0.f);
#pragma omp section
		restv.init(m_restv.getXSize(), m_restv.getYSize(), m_restv.getZSize(), 0.f);
#pragma omp section
		restw.init(m_restw.getXSize(), m_restw.getYSize(), m_restw.getZSize(), 0.f);
	}

	if (m_enableFire && m_advectFuel && m_fuelAdvectRate)
	    fuel.init(m_fuel.getXSize(), m_fuel.getYSize(), m_fuel.getZSize(), 0.f);

	if (m_enableFire)
		heat.init(m_heat.getXSize(), m_heat.getYSize(), m_heat.getZSize(), 0.f);

	// update smoke and temperature 
#pragma omp parallel for    
	for (int i = 0; i < m_pressure.getXSize(); i++)
	{
		for (int j = 0; j < m_pressure.getYSize(); j++)
		{
			for (int k = 0; k < m_pressure.getZSize(); k++)
			{
				if (m_layer(i, j, k) != SOLID && m_layer(i, j, k) == DENSITY)
				{
					float x, y, z;
					vector3 newpnt;
					float nusmoke, nutemp, nufuel, nuheat;
					vector3 nurest;

					x = i * m_dx;
					y = j * m_dx;
					z = k * m_dx;

					switch (m_advectMode)
					{
					case BFECC:
						newpnt = bfecc(x, y, z, -m_timeStep, 0.5f);
						break;
					default:
						newpnt = traceParticle(x, y, z, -m_timeStep, 0.5f);
					}
					
					nusmoke = getSmoke(newpnt.x, newpnt.y, newpnt.z);
					nutemp = getTemperature(newpnt.x, newpnt.y, newpnt.z);
					nurest = getRest(newpnt.x, newpnt.y, newpnt.z);
					smoke.set(nusmoke, i, j, k);
					temperature.set(nutemp, i, j, k);
					restu.set(nurest.x, i, j, k);
					restv.set(nurest.y, i, j, k);
					restw.set(nurest.z, i, j, k);

					if (m_enableFire)
					{
						nuheat = getHeat(newpnt.x, newpnt.y, newpnt.z);
						heat.set(nuheat, i, j, k);
					}

					// fuel needs to be advected separately
					if (m_enableFire && m_advectFuel && m_fuelAdvectRate)
					{
						switch (m_advectMode)
						{
						case BFECC:
							newpnt = bfecc(x, y, z, -m_timeStep * m_fuelAdvectRate, 0.f);
							break;
						default:
							newpnt = traceParticle(x, y, z, -m_timeStep * m_fuelAdvectRate, 0.f);
						}

						nufuel = getFuel(newpnt.x, newpnt.y, newpnt.z);
						fuel.set(nufuel, i, j, k);
					}
				}
			}
		}
	}

#pragma omp parallel for    
	for (int i = 0; i < m_smoke.getNumVoxels(); i++)
	{
		m_smoke.setIndexValue(i, smoke[i]);
		m_temperature.setIndexValue(i, temperature[i]);
		m_restu.setIndexValue(i, restu[i]);
		m_restv.setIndexValue(i, restv[i]);
		m_restw.setIndexValue(i, restw[i]);

		if (m_enableFire)
		    m_heat.setIndexValue(i, heat[i]);

		if (m_advectFuel && m_fuelAdvectRate)
		{
			m_fuel.setIndexValue(i, fuel[i]);
		}
	}
}

void GasSolver::applyBouyancyForce()
{
	float alpha = m_smokeWeight;
	float beta = m_bouyancyLift;

	float step = 0.5 * m_timeStep;
#pragma omp parallel for    
	for (int i = 0; i < m_pressure.getXSize(); i++)
	{
		for (int j = 0; j < m_pressure.getYSize(); j++)
		{
			for (int k = 0; k < m_pressure.getZSize(); k++)
			{
				if (m_layer(i, j, k) != SOLID)
				{
					float velu, velv, velw;
					float smoke = m_smoke(i, j, k);
					float temperature = m_temperature(i, j, k);
					vector3 bouyancy;
					bouyancy.x = step * (alpha * smoke - beta * (temperature - m_ambientTemperature)) * m_gravity.x;
					bouyancy.y = step * (alpha * smoke - beta * (temperature - m_ambientTemperature)) * m_gravity.y;
					bouyancy.z = step * (alpha * smoke - beta * (temperature - m_ambientTemperature)) * m_gravity.z;

					velu = m_velu(i, j, k) + bouyancy.x;
					m_velu.set(velu, i, j, k);
					velu = m_velu(i + 1, j, k) + bouyancy.x;
					m_velu.set(velu, i + 1, j, k);

					velv = m_velv(i, j, k) + bouyancy.y;
					m_velv.set(velv, i, j, k);
					velv = m_velv(i, j + 1, k) + bouyancy.y;
					m_velv.set(velv, i, j + 1, k);

					velw = m_velw(i, j, k) + bouyancy.z;
					m_velw.set(velw, i, j, k);
					velw = m_velw(i, j, k + 1) + bouyancy.z;
					m_velw.set(velw, i, j, k + 1);
				}
			}
		}
	}
}

// this implementation is unstable as (m_kinematicViscosity * m_timeStep > 1/6)
// need to solve iteratively 
void GasSolver::applyViscosity()
{
	Grid3D<float> velu;
	Grid3D<float> velv;
	Grid3D<float> velw;

#pragma omp parallel sections
	{
#pragma omp section
		velu.init(m_velu.getXSize(), m_velu.getYSize(), m_velu.getZSize(), 0.f);
#pragma omp section
		velv.init(m_velv.getXSize(), m_velv.getYSize(), m_velv.getZSize(), 0.f);
#pragma omp section
		velw.init(m_velw.getXSize(), m_velw.getYSize(), m_velw.getZSize(), 0.f);
	}

	// update u
#pragma omp parallel for   
	for (int i = 0; i < m_velu.getXSize(); i++)
	{
		for (int j = 0; j < m_velu.getYSize(); j++)
		{
			for (int k = 0; k < m_velu.getZSize(); k++)
			{
				float vu = m_velu(i + 1, j, k) + m_velu(i - 1, j, k) + m_velu(i, j + 1, k) + m_velu(i, j - 1, k) + m_velu(i, j, k + 1) + m_velu(i, j, k - 1) - 6 * m_velu(i, j, k);
				velu.set(vu, i, j, k);
			}
		}
	}

#pragma omp parallel for    
	for (int i = 0; i < m_velu.getNumVoxels(); i++)
	{
		float newvel = m_velu[i] + m_kinematicViscosity * m_timeStep * velu[i];
		m_velu.setIndexValue(i, newvel);
	}

	// update v
#pragma omp parallel for    
	for (int i = 0; i < m_velv.getXSize(); i++)
	{
		for (int j = 0; j < m_velv.getYSize(); j++)
		{
			for (int k = 0; k < m_velv.getZSize(); k++)
			{
				float vv = m_velv(i + 1, j, k) + m_velv(i - 1, j, k) + m_velv(i, j + 1, k) + m_velv(i, j - 1, k) + m_velv(i, j, k + 1) + m_velv(i, j, k - 1) - 6 * m_velv(i, j, k);
				velv.set(vv, i, j, k);
			}
		}
	}

#pragma omp parallel for    
	for (int i = 0; i < m_velv.getNumVoxels(); i++)
	{
		float newvel = m_velv[i] + m_kinematicViscosity * m_timeStep * velv[i];
		m_velv.setIndexValue(i, newvel);
	}

	// update w
#pragma omp parallel for    
	for (int i = 0; i < m_velw.getXSize(); i++)
	{
		for (int j = 0; j < m_velw.getYSize(); j++)
		{
			for (int k = 0; k < m_velw.getZSize(); k++)
			{
				float vw = m_velw(i + 1, j, k) + m_velw(i - 1, j, k) + m_velw(i, j + 1, k) + m_velw(i, j - 1, k) + m_velw(i, j, k + 1) + m_velw(i, j, k - 1) - 6 * m_velw(i, j, k);
				velw.set(vw, i, j, k);
			}
		}
	}

#pragma omp parallel for    
	for (int i = 0; i < m_velw.getNumVoxels(); i++)
	{
		float newvel = m_velw[i] + m_kinematicViscosity * m_timeStep * velw[i];
		m_velw.setIndexValue(i, newvel);
	}
}

void GasSolver::project()
{
	// Solve  A * p = b

	// construct A matrix 
#pragma omp parallel sections
	{
#pragma omp section
		if (m_computePressureCoefficients)
		{
			computePressureCoefficients();
		}

		// construct B matrix  
#pragma omp section
		computeNegativeDivergence();
	}

	// burn fuel (this modifies divergence so doing this step here?)
	if (m_enableFire)
		burnFuel();

	// solve for p
	//PressureSolver solver;
	//solver.setMaxIterations(m_maxPressureSolveIterations);
	//solver.setTolerance(m_pressureSolveTolerance);
	int nvoxels = m_pressure.getNumVoxels();
	VectorXd x;
	x.resize(nvoxels);
	x.setZero();
	
	if (m_pressureSolver == MULTIGRID)
	{
		//x = solver.jacobiSolver(true, m_coefficients, m_divergence, m_guess);
		solver.initMultiGrid(5, m_pressure.getXSize(), m_pressure.getYSize(), m_pressure.getZSize());
		x = solver.multiGridSolver_Vcycle(true, m_coefficients, m_divergence, m_guess);
	}
	else
	{
		x = solver.conjugateGradientSolver(m_coefficients, m_divergence, m_guess);
	}

#pragma omp parallel for   
	for (int i = 0; i < nvoxels; i++)
	{
		m_guess[i] = x[i];
	}

	m_pressure.setZero();

#pragma omp parallel for    
	for (int i = 0; i < nvoxels; i++)
	{
		m_pressure.setIndexValue(i, x(i));
	}

	// apply pressure
	pressureUpdate();
}

void GasSolver::computePressureCoefficients()
{
	int nvoxels = m_pressure.getNumVoxels();
	float rho = 1.0;
	float scale = m_timeStep / (rho * m_dx * m_dx);

	m_coefficients.setZero();
	m_coefficients.reserve(VectorXi::Constant(nvoxels, 7));

#pragma omp parallel for    
	for (int i = 0; i < m_pressure.getXSize(); i++)
	{
		for (int j = 0; j < m_pressure.getYSize(); j++)
		{
			for (int k = 0; k < m_pressure.getZSize(); k++)
			{
				int IX = i * m_pressure.getYSize() *  m_pressure.getZSize() + j *  m_pressure.getZSize() + k;
				int IXPI = (i + 1) * m_pressure.getYSize() *  m_pressure.getZSize() + j *  m_pressure.getZSize() + k;
				int IXPJ = i * m_pressure.getYSize() *  m_pressure.getZSize() + (j + 1) *  m_pressure.getZSize() + k;
				int IXPK = i * m_pressure.getYSize() *  m_pressure.getZSize() + j *  m_pressure.getZSize() + (k + 1);
				int IXNI = (i - 1) * m_pressure.getYSize() *  m_pressure.getZSize() + j *  m_pressure.getZSize() + k;
				int IXNJ = i * m_pressure.getYSize() *  m_pressure.getZSize() + (j - 1) *  m_pressure.getZSize() + k;
				int IXNK = i * m_pressure.getYSize() *  m_pressure.getZSize() + j *  m_pressure.getZSize() + (k - 1);

				if (m_layer(i, j, k) != SOLID && m_layer(i + 1, j, k) != SOLID)
				{
					m_coefficients.coeffRef(IX, IX) += scale;
					if (IXPI < nvoxels)
					{
						m_coefficients.coeffRef(IX, IXPI) = -scale;
					}
				}

				if (m_layer(i, j, k) != SOLID && m_layer(i - 1, j, k) != SOLID)
				{
					m_coefficients.coeffRef(IX, IX) += scale;
					if (IXNI >= 0)
					{
						m_coefficients.coeffRef(IX, IXNI) = -scale;
					}
				}

				if (m_layer(i, j, k) != SOLID && m_layer(i, j + 1, k) != SOLID)
				{
					m_coefficients.coeffRef(IX, IX) += scale;
					if (IXPJ < nvoxels)
					{
					    m_coefficients.coeffRef(IX, IXPJ) = -scale;
					}
				}

				if (m_layer(i, j, k) != SOLID && m_layer(i, j - 1, k) != SOLID)
				{
					m_coefficients.coeffRef(IX, IX) += scale;
					if (IXNJ >= 0)
					{
						m_coefficients.coeffRef(IX, IXNJ) = -scale;
					}
				}

				if (m_layer(i, j, k) != SOLID && m_layer(i, j, k + 1) != SOLID)
				{
					m_coefficients.coeffRef(IX, IX) += scale;
					if (IXPK < nvoxels)
					{
						m_coefficients.coeffRef(IX, IXPK) = -scale;
					}
				}

				if (m_layer(i, j, k) != SOLID && m_layer(i, j, k - 1) != SOLID)
				{
					m_coefficients.coeffRef(IX, IX) += scale;
					if (IXNK >= 0)
					{
						m_coefficients.coeffRef(IX, IXNK) = -scale;
					}
				}
			}
		}
	}

	m_coefficients.makeCompressed();
}

void GasSolver::computeNegativeDivergence()
{
	float scale = 1.0 / m_dx;   

	m_divergence.setZero();
#pragma omp parallel for    
	for (int i = 0; i < m_pressure.getXSize(); i++)
	{
		for (int j = 0; j < m_pressure.getYSize(); j++)
		{
			for (int k = 0; k < m_pressure.getZSize(); k++)
			{
				if (m_layer(i, j, k) != SOLID)
				{
					float val = -scale * divergence(i, j, k);
					int IX = i * m_pressure.getYSize() *  m_pressure.getZSize() + j *  m_pressure.getZSize() + k;
					m_divergence(IX) = val;
				}
			}
		}
	}

	// account for solid velocities
#pragma omp parallel for    
	for (int i = 0; i < m_pressure.getXSize(); i++)
	{
		for (int j = 0; j < m_pressure.getYSize(); j++)
		{
			for (int k = 0; k < m_pressure.getZSize(); k++)
			{
				int IX = i * m_pressure.getYSize() *  m_pressure.getZSize() + j *  m_pressure.getZSize() + k;

				if (m_layer(i, j, k) != SOLID)
				{
					if (m_layer(i - 1, j, k) == SOLID)
					{
						float val = m_divergence(IX) - scale * (m_velu(i, j, k) - m_solidu(i, j, k));
						m_divergence(IX) = val;
					}
					if (m_layer(i + 1, j, k) == SOLID)
					{
						float val = m_divergence(IX) + scale * (m_velu(i + 1, j, k) - m_solidu(i + 1, j, k));
						m_divergence(IX) = val;
					}

					if (m_layer(i, j - 1, k) == SOLID)
					{
						float val = m_divergence(IX) - scale * (m_velv(i, j, k) - m_solidv(i, j, k));
						m_divergence(IX) = val;
					}
					if (m_layer(i, j + 1, k) == SOLID)
					{
						float val = m_divergence(IX) + scale * (m_velv(i, j + 1, k) - m_solidv(i, j + 1, k));
						m_divergence(IX) = val;
					}

					if (m_layer(i, j, k - 1) == SOLID)
					{
						float val = m_divergence(IX) - scale * (m_velw(i, j, k) - m_solidw(i, j, k));
						m_divergence(IX) = val;
					}
					if (m_layer(i, j, k + 1) == SOLID)
					{
						float val = m_divergence(IX) + scale * (m_velw(i, j, k + 1) - m_solidw(i, j, k + 1));
						m_divergence(IX) = val;
					}
				}
			}
		}
	}
}

void GasSolver::pressureUpdate()
{
	float newvel;
	float rho = 1.0;
	float scale = m_timeStep / (rho * m_dx);

#pragma omp parallel for    
	for (int i = 0; i < m_pressure.getXSize(); i++)
	{
		for (int j = 0; j < m_pressure.getYSize(); j++)
		{
			for (int k = 0; k < m_pressure.getZSize(); k++)
			{
				if (m_layer(i, j, k) != SOLID)
				{
					newvel = m_velu(i, j, k) - scale * m_pressure(i, j, k);
					m_velu.set(newvel, i, j, k);
					newvel = m_velu(i + 1, j, k) + scale * m_pressure(i, j, k);
					m_velu.set(newvel, i + 1, j, k);

					newvel = m_velv(i, j, k) - scale * m_pressure(i, j, k);
					m_velv.set(newvel, i, j, k);
					newvel = m_velv(i, j + 1, k) + scale * m_pressure(i, j, k);
					m_velv.set(newvel, i, j + 1, k);

					newvel = m_velw(i, j, k) - scale * m_pressure(i, j, k);
					m_velw.set(newvel, i, j, k);
					newvel = m_velw(i, j, k + 1) + scale * m_pressure(i, j, k);
					m_velw.set(newvel, i, j, k + 1);
				}
			}
		}
	}

#pragma omp parallel for    
	for (int i = 0; i < m_pressure.getXSize(); i++)
	{
		for (int j = 0; j < m_pressure.getYSize(); j++)
		{
			for (int k = 0; k < m_pressure.getZSize(); k++)
			{
			    if (m_layer(i, j, k) == SOLID)
				{
					newvel = m_solidu(i, j, k);
					m_velu.set(newvel, i, j, k);
					newvel = m_solidu(i + 1, j, k);
					m_velu.set(newvel, i + 1, j, k);

					newvel = m_solidv(i, j, k);
					m_velv.set(newvel, i, j, k);
					newvel = m_solidv(i, j + 1, k);
					m_velv.set(newvel, i, j + 1, k);

					newvel = m_solidw(i, j, k);
					m_velw.set(newvel, i, j, k);
					newvel = m_solidw(i, j, k + 1);
					m_velw.set(newvel, i, j, k + 1);
				}
			}
		}
	}
}

float GasSolver::divergence(int ix, int iy, int iz)
{
	float retval;
    retval = m_velu(ix + 1, iy, iz) - m_velu(ix, iy, iz) + m_velv(ix, iy + 1, iz) - m_velv(ix, iy, iz) + m_velw(ix, iy, iz + 1) - m_velw(ix, iy, iz);
	return retval;
}

vector3 GasSolver::pressureGradient(int i, int j, int k)
{
	vector3 retval;

	retval.x = m_pressure(i + 1, j, k) - m_pressure(i, j, k);
	retval.y = m_pressure(i, j + 1, k) - m_pressure(i, j, k);
	retval.z = m_pressure(i, j, k + 1) - m_pressure(i, j, k);

	return retval;
}

vector3 GasSolver::temperatureGradient(int i, int j, int k)
{
	vector3 retval;

	retval.x = m_temperature(i + 1, j, k) - m_temperature(i, j, k);
	retval.y = m_temperature(i, j + 1, k) - m_temperature(i, j, k);
	retval.z = m_temperature(i, j, k + 1) - m_temperature(i, j, k);

	return retval;
}

vector3 GasSolver::heatGradient(int i, int j, int k)
{
	vector3 retval;

	retval.x = m_heat(i + 1, j, k) - m_heat(i, j, k);
	retval.y = m_heat(i, j + 1, k) - m_heat(i, j, k);
	retval.z = m_heat(i, j, k + 1) - m_heat(i, j, k);

	return retval;
}

void GasSolver::applyBoundaryConditions()
{
#pragma omp parallel for    
	for (int j = 0; j < m_velu.getYSize(); j++)
	{
		for (int k = 0; k < m_velu.getZSize(); k++)
		{
			float velu = m_velu(3, j, k);
			m_velu.set(velu, 2, j, k);
			m_velu.set(velu, 1, j, k);
			m_velu.set(velu, 0, j, k);
			velu = m_velu(m_velu.getXSize() - 4, j, k);
			m_velu.set(velu, m_velu.getXSize() - 3, j, k);
			m_velu.set(velu, m_velu.getXSize() - 2, j, k);
			m_velu.set(velu, m_velu.getXSize() - 1, j, k);

			float velv = m_velv(2, j, k);
			m_velv.set(velv, 1, j, k);
			m_velv.set(velv, 0, j, k);			
			velv = m_velv(m_velv.getXSize() - 3, j, k);
			m_velv.set(velv, m_velv.getXSize() - 2, j, k);
			m_velv.set(velv, m_velv.getXSize() - 1, j, k);

			float velw = m_velw(2, j, k);
			m_velw.set(velw, 1, j, k);
			m_velw.set(velw, 0, j, k);
			velw = m_velw(m_velw.getXSize() - 3, j, k);
			m_velw.set(velw, m_velw.getXSize() - 2, j, k);
			m_velw.set(velw, m_velw.getXSize() - 1, j, k);
		}
	}

#pragma omp parallel for    
	for (int i = 0; i < m_velv.getXSize(); i++)
	{
		for (int k = 0; k < m_velv.getZSize(); k++)
		{
			float velv = m_velv(i, 3, k);
			m_velv.set(velv, i, 2, k);
			m_velv.set(velv, i, 1, k);
			m_velv.set(velv, i, 0, k);
			velv = m_velv(i, m_velv.getYSize() - 4, k);
			m_velv.set(velv, i, m_velv.getYSize() - 3, k);
			m_velv.set(velv, i, m_velv.getYSize() - 2, k);
			m_velv.set(velv, i, m_velv.getYSize() - 1, k);

			float velu = m_velu(i, 2, k);
			m_velu.set(velu, i, 1, k);
			m_velu.set(velu, i, 0, k);
			velu = m_velu(i, m_velu.getYSize() - 3, k);
			m_velu.set(velu, i, m_velu.getYSize() - 2, k);
			m_velu.set(velu, i, m_velu.getYSize() - 1, k);

			float velw = m_velw(i, 2, k);
			m_velw.set(velw, i, 1, k);
			m_velw.set(velw, i, 0, k);
			velw = m_velw(i, m_velw.getYSize() - 3, k);
			m_velw.set(velw, i, m_velw.getYSize() - 2, k);
			m_velw.set(velw, i, m_velw.getYSize() - 1, k);
		}
	}

#pragma omp parallel for    
	for (int i = 0; i < m_velw.getXSize(); i++)
	{
		for (int j = 0; j < m_velw.getYSize(); j++)
		{
			float velw = m_velw(i, j, 3);
			m_velw.set(velw, i, j, 2);
			m_velw.set(velw, i, j, 1);
			m_velw.set(velw, i, j, 0.f);
			velw = m_velw(i, j, m_velw.getZSize() - 4);
			m_velw.set(velw, i, j, m_velw.getZSize() - 3);
			m_velw.set(velw, i, j, m_velw.getZSize() - 2);
			m_velw.set(velw, i, j, m_velw.getZSize() - 1);

			float velu = m_velu(i, j, 2);
			m_velu.set(velu, i, j, 1);
			m_velu.set(velu, i, j, 0.f);
			velu = m_velu(i, j, m_velu.getZSize() - 3);
			m_velu.set(velu, i, j, m_velu.getZSize() - 2);
			m_velu.set(velu, i, j, m_velu.getZSize() - 1);

			float velv = m_velv(i, j, 2);
			m_velv.set(velv, i, j, 1);
			m_velv.set(velv, i, j, 0.f);
			velv = m_velv(i, j, m_velv.getZSize() - 3);
			m_velv.set(velv, i, j, m_velv.getZSize() - 2);
			m_velv.set(velv, i, j, m_velv.getZSize() - 1);
		}
	}
}

// velocity at cell center
vector3 GasSolver::cellCenterVelocity(int ix, int iy, int iz)
{
	vector3 vavg;

	vavg.x = 0.5f * (m_velu(ix, iy, iz) + m_velu(ix + 1, iy, iz));
	vavg.y = 0.5f * (m_velv(ix, iy, iz) + m_velv(ix, iy + 1, iz));
	vavg.z = 0.5f * (m_velw(ix, iy, iz) + m_velw(ix, iy, iz + 1));

	return vavg;
}

// velocities at cell faces
vector3 GasSolver::faceCenterVelocity(int face, int ix, int iy, int iz)
{
	vector3 vavg = { 0,0,0 };
	vector3 velxp, velyp, velzp;

	velxp.x = m_velu(ix + 1, iy, iz);
	velxp.y = 0.25f * (m_velv(ix, iy, iz) + m_velv(ix, iy + 1, iz) + m_velv(ix + 1, iy, iz) + m_velv(ix + 1, iy + 1, iz));
	velxp.z = 0.25f * (m_velw(ix, iy, iz) + m_velw(ix, iy, iz + 1) + m_velw(ix + 1, iy, iz) + m_velw(ix + 1, iy, iz + 1));

	velyp.x = 0.25f * (m_velu(ix, iy, iz) + m_velu(ix + 1, iy, iz) + m_velu(ix, iy + 1, iz) + m_velu(ix + 1, iy + 1, iz));
	velyp.y = m_velv(ix, iy + 1, iz);
	velyp.z = 0.25f * (m_velw(ix, iy, iz) + m_velw(ix, iy, iz + 1) + m_velw(ix, iy + 1, iz) + m_velw(ix, iy + 1, iz + 1));

	velzp.x = 0.25f * (m_velu(ix, iy, iz) + m_velu(ix + 1, iy, iz) + m_velu(ix, iy, iz + 1) + m_velu(ix + 1, iy, iz + 1));
	velzp.y = 0.25f * (m_velv(ix, iy, iz) + m_velv(ix, iy + 1, iz) + m_velv(ix, iy, iz + 1) + m_velv(ix, iy + 1, iz + 1));
	velzp.z = m_velw(ix, iy, iz + 1);

	switch (face)
	{ 
	case 0:
		vavg.x = velxp.x; vavg.y = velxp.y; vavg.z = velxp.z;
		return vavg;
	case 1:
		vavg.x = velyp.x; vavg.y = velyp.y; vavg.z = velyp.z;
		return vavg;
	case 2:
		vavg.x = velzp.x; vavg.y = velzp.y; vavg.z = velzp.z;
		return vavg;
	}

	return vavg;
}

void GasSolver::applyVorticityConfinement()
{
	Grid3D<vector3> curlVec;
	Grid3D<vector3> gradCurl;

	vector3 zerovec = { 0.f, 0.f, 0.f };

#pragma omp parallel sections
	{
#pragma omp section
		curlVec.init(m_xsize, m_ysize, m_zsize, zerovec);
#pragma omp section
		gradCurl.init(m_xsize, m_ysize, m_zsize, zerovec);
	}

	// compute curl
#pragma omp parallel for    
	for (int i = 0; i < m_pressure.getXSize(); i++)
	{
		for (int j = 0; j < m_pressure.getYSize(); j++)
		{
			for (int k = 0; k < m_pressure.getZSize(); k++)
			{
				if (m_layer(i, j, k) != SOLID)
				{
					vector3 curl;
					curl.x = (cellCenterVelocity(i, j + 1, k).z - cellCenterVelocity(i, j - 1, k).z) / (2.0 * m_dx) -
						(cellCenterVelocity(i, j, k + 1).y - cellCenterVelocity(i, j, k - 1).y) / (2.0 * m_dx);
					curl.y = (cellCenterVelocity(i, j, k + 1).x - cellCenterVelocity(i, j, k - 1).x) / (2.0 * m_dx) -
						(cellCenterVelocity(i + 1, j, k).z - cellCenterVelocity(i - 1, j, k).z) / (2.0 * m_dx);
					curl.z = (cellCenterVelocity(i + 1, j, k).y - cellCenterVelocity(i - 1, j, k).y) / (2.0 * m_dx) -
						(cellCenterVelocity(i, j + 1, k).x - cellCenterVelocity(i, j - 1, k).x) / (2.0 * m_dx);

					curlVec.set(curl, i, j, k);
				}
			}
		}
	}

	// compute gradient of curl
#pragma omp parallel for    
	for (int i = 0; i < m_pressure.getXSize(); i++)
	{
		for (int j = 0; j < m_pressure.getYSize(); j++)
		{
			for (int k = 0; k < m_pressure.getZSize(); k++)
			{
				if (m_layer(i, j, k) != SOLID)
				{
					float curlmagxp, curlmagxn, curlmagyp, curlmagyn, curlmagzp, curlmagzn;

					curlmagxp = sqrt(curlVec(i + 1, j, k).x * curlVec(i + 1, j, k).x + curlVec(i + 1, j, k).y * curlVec(i + 1, j, k).y + 
						curlVec(i + 1, j, k).z * curlVec(i + 1, j, k).z);
					curlmagxn = sqrt(curlVec(i - 1, j, k).x * curlVec(i - 1, j, k).x + curlVec(i - 1, j, k).y * curlVec(i - 1, j, k).y + 
						curlVec(i - 1, j, k).z * curlVec(i - 1, j, k).z);

					curlmagyp = sqrt(curlVec(i, j + 1, k).x * curlVec(i, j + 1, k).x + curlVec(i, j + 1, k).y * curlVec(i, j + 1, k).y + 
						curlVec(i, j + 1, k).z * curlVec(i, j + 1, k).z);
					curlmagyn = sqrt(curlVec(i, j - 1, k).x * curlVec(i, j - 1, k).x + curlVec(i, j - 1, k).y * curlVec(i, j - 1, k).y + 
						curlVec(i, j - 1, k).z * curlVec(i, j - 1, k).z);

					curlmagzp = sqrt(curlVec(i, j, k + 1).x * curlVec(i, j, k + 1).x + curlVec(i, j, k + 1).y * curlVec(i, j, k + 1).y + 
						curlVec(i, j, k + 1).z * curlVec(i, j, k + 1).z);
					curlmagzn = sqrt(curlVec(i, j, k - 1).x * curlVec(i, j, k - 1).x + curlVec(i, j, k - 1).y * curlVec(i, j, k - 1).y + 
						curlVec(i, j, k - 1).z * curlVec(i, j, k - 1).z);

					vector3 gradc;
					gradc.x = (curlmagxp - curlmagxn) / (2.0 * m_dx);
					gradc.y = (curlmagyp - curlmagyn) / (2.0 * m_dx);
					gradc.z = (curlmagzp - curlmagzn) / (2.0 * m_dx);

					float gradmag;
					float zerodiv = 1.0e-20;
					gradmag = sqrt(gradc.x * gradc.x + gradc.y * gradc.y + gradc.z * gradc.z) + (zerodiv / (m_dx * m_timeStep));
					gradc.x /= gradmag;
					gradc.y /= gradmag;
					gradc.z /= gradmag;

					gradCurl.set(gradc, i, j, k);
				}
			}
		}
	}

	// compute vorticity confinement force 
	float step = 0.5f * m_timeStep * m_vorticityScale * m_dx;
#pragma omp parallel for    
	for (int i = 0; i < m_pressure.getXSize(); i++)
	{
		for (int j = 0; j < m_pressure.getYSize(); j++)
		{
			for (int k = 0; k < m_pressure.getZSize(); k++)
			{
				if (m_layer(i, j, k) != SOLID)
				{

					vector3 xprod;
					xprod.x = gradCurl(i, j, k).z * curlVec(i, j, k).y - gradCurl(i, j, k).y * curlVec(i, j, k).z;
					xprod.y = gradCurl(i, j, k).x * curlVec(i, j, k).z - gradCurl(i, j, k).z * curlVec(i, j, k).x;
					xprod.z = gradCurl(i, j, k).y * curlVec(i, j, k).x - gradCurl(i, j, k).x * curlVec(i, j, k).y;

					float velu = m_velu(i, j, k) + step * xprod.x;
					m_velu.set(velu, i, j, k);

					velu = m_velu(i + 1, j, k) + step * xprod.x;
					m_velu.set(velu, i + 1, j, k);

					float velv = m_velv(i, j, k) + step * xprod.y;
					m_velv.set(velv, i, j, k);

					velv = m_velv(i, j + 1, k) + step * xprod.y;
					m_velv.set(velv, i, j + 1, k);

					float velw = m_velw(i, j, k) + step * xprod.z;
					m_velw.set(velw, i, j, k);
				
					velw = m_velw(i, j, k + 1) + step * xprod.z;
					m_velw.set(velw, i, j, k + 1);
				}
			}
		}
	}
}

void GasSolver::applyWindForce()
{
	int nvoxels;
	float velu, velv, velw;

	nvoxels = m_velu.getNumVoxels();
#pragma omp parallel for    
	for (int i = 0; i < nvoxels; i++) 
	{
		velu = (1.0f - m_kwind * m_timeStep) * m_velu[i] + m_kwind * m_timeStep * m_windForce.x;
		m_velu.setIndexValue(i, velu);
	}

	nvoxels = m_velv.getNumVoxels();
#pragma omp parallel for    
	for (int i = 0; i < nvoxels; i++)
	{
		velv = (1.0f - m_kwind * m_timeStep) * m_velv[i] + m_kwind * m_timeStep * m_windForce.y;
		m_velv.setIndexValue(i, velv);
	}

	nvoxels = m_velw.getNumVoxels();
#pragma omp parallel for    
	for (int i = 0; i < nvoxels; i++)
	{
		velw = (1.0f - m_kwind * m_timeStep) * m_velw[i] + m_kwind * m_timeStep * m_windForce.z;
		m_velw.setIndexValue(i, velw);
	}
}

void GasSolver::applyDiffusion()
{
#pragma omp parallel for    
	for (int i = 0; i < m_pressure.getXSize(); i++)
	{
		for (int j = 0; j < m_pressure.getYSize(); j++)
		{
			for (int k = 0; k < m_pressure.getZSize(); k++)
			{
				// Gaussian blur (m_diffusion = radius of kernel), multiply by m_dt
				for (int l = 0; l < m_diffusion; k++)
				{
				}
			}
		}
	}
}

void GasSolver::applyDissipation()
{
#pragma omp parallel for    
	for (int i = 0; i < m_pressure.getXSize(); i++)
	{
		for (int j = 0; j < m_pressure.getYSize(); j++)
		{
			for (int k = 0; k < m_pressure.getZSize(); k++)
			{
				if (m_layer(i, j, k) != SOLID)
				{
					float decay = exp(-m_dissipation * m_timeStep) * m_smoke(i, j, k);
					m_smoke.set(decay, i, j, k);
				}
			}
		}
	}
}

void GasSolver::applyCooling()
{
#pragma omp parallel for    
	for (int i = 0; i < m_pressure.getXSize(); i++)
	{
		for (int j = 0; j < m_pressure.getYSize(); j++)
		{
			for (int k = 0; k < m_pressure.getZSize(); k++)
			{
				if (m_layer(i, j, k) != SOLID)
				{
					float decay = exp(-m_cooling * m_timeStep) * m_temperature(i, j, k);
					m_temperature.set(decay, i, j, k);
				}
			}
		}
	}
}

void GasSolver::updateFluidSources()
{
	for (std::vector<ImplicitSphereSource>::iterator it = m_implicitSphereSources.begin(); it != m_implicitSphereSources.end(); ++it)
	{
		if (it->isActive())
		{
			it->update(m_smoke, m_temperature, m_fuel, m_velu, m_velv, m_velw);
		}
	}

	for (std::vector<GeometrySource>::iterator it = m_geometrySources.begin(); it != m_geometrySources.end(); ++it)
	{
		if (it->isActive())
		{
			it->update(m_smoke, m_temperature, m_fuel, m_velu, m_velv, m_velw);
		}
	}
}

void GasSolver::updateColliders()
{
	bool updatePressureCoefficients = false;

	for (std::vector<ImplicitSphereCollider>::iterator it = m_implicitSphereColliders.begin(); it != m_implicitSphereColliders.end(); ++it)
	{
		if (it->isActive())
		{
			if (!it->getStaticFlag())
			{
				updatePressureCoefficients = true;
			}
		}
	}

	for (std::vector<GeometryCollider>::iterator it = m_geometryColliders.begin(); it != m_geometryColliders.end(); ++it)
	{
		if (it->isActive())
		{
			if (!it->getStaticFlag())
			{
				updatePressureCoefficients = true;
			}
		}
	}

	m_computePressureCoefficients = updatePressureCoefficients;
}

void GasSolver::applyTurbulence()
{
	float maxdensity = 1.0f;

	vector3 center = { 0.5f * (m_pressure.getXSize() - 0.5f) * m_dx, 0.5f * (m_pressure.getYSize() - 0.5f) * m_dx, 0.5f * (m_pressure.getZSize() - 0.5f) * m_dx };

	float step = 0.5f * m_timeStep * m_fps * m_turbulence;

#pragma omp parallel for    
	for (int i = 0; i < m_pressure.getXSize(); i++)
	{
		for (int j = 0; j < m_pressure.getYSize(); j++)
		{
			for (int k = 0; k < m_pressure.getZSize(); k++)
			{
				if (m_layer(i, j, k) != SOLID)
				{
					float density = m_smoke(i, j, k);
					float heat = m_heat(i, j, k);
					if (density > m_turbulenceThreshold || heat > m_turbulenceThreshold)
					{
						vector3 rest;
						float x = i * m_dx;
						float y = j * m_dx;
						float z = k * m_dx;
						rest = getRest(x, y, z);
						vector3 turb = m_curlNoise.getCurlVelocity(rest, center, m_flow * m_simTime);
						turb = turb * step;

						float controlField = 1.0f;
						if (m_enableTurbulenceControlField)
							controlField *= 1.0f - smoothstep(m_turbulenceThreshold, maxdensity, density);

						float velu = m_velu(i, j, k) + controlField * turb.x;
					    m_velu.set(velu, i, j, k);

						velu = m_velu(i + 1, j, k) + controlField * turb.x;
						m_velu.set(velu, i + 1, j, k);

						float velv = m_velv(i, j, k) + controlField * turb.y;
					    m_velv.set(velv, i, j, k);

						velv = m_velv(i, j + 1, k) + controlField * turb.y;
						m_velv.set(velv, i, j + 1, k);

						float velw = m_velw(i, j, k) + controlField * turb.z;
					    m_velw.set(velw, i, j, k);

						velw = m_velw(i, j, k + 1) + controlField * turb.z;
						m_velw.set(velw, i, j, k + 1);
					}
				}
			}
		}
	}
}

void GasSolver::applyFuelDissipation()
{
#pragma omp parallel for    
	for (int i = 0; i < m_pressure.getXSize(); i++)
	{
		for (int j = 0; j < m_pressure.getYSize(); j++)
		{
			for (int k = 0; k < m_pressure.getZSize(); k++)
			{
				if (m_layer(i, j, k) != SOLID)
				{
					float decay = exp(-m_fuelDissipation * m_timeStep) * m_fuel(i, j, k);
					if (decay < 0) decay = 0;
					m_fuel.set(decay, i, j, k);
				}
			}
		}
	}
}

void GasSolver::burnFuel()
{
#pragma omp parallel for    
	for (int i = 0; i < m_pressure.getXSize(); i++)
	{
		for (int j = 0; j < m_pressure.getYSize(); j++)
		{
			for (int k = 0; k < m_pressure.getZSize(); k++)
			{
				if (m_layer(i, j, k) != SOLID)
				{
					if (m_temperature(i, j, k) >= m_ignitionTemperature && m_fuel(i, j, k) > 0)
					{
						float fuel = m_fuel(i, j, k);
						float nufuel = m_fuel(i, j, k) - m_burnRate * m_timeStep;
						if (nufuel < 0) nufuel = 0;
						m_fuel.set(nufuel, i, j, k);
						float deltaFuel = fuel - nufuel;

						// need better temperature model, should be both an ignition AND max temperature? see book pg. 105.
						// update temperature and smoke based on burnt fuel (heat)

						float heat = m_heat(i, j, k) + m_fuelBurnRate * deltaFuel;
						m_heat.set(heat, i, j, k);

						if (m_enableFireSmoke && heat > m_heatCutoff)
						{
							if (m_enableDenseSmoke)
							{
								float smoke = m_smoke(i, j, k) + m_fuelSmokeRate * m_heat(i, j, k);
								m_smoke.set(smoke, i, j, k);
							}
							else
							{
								float smoke = m_smoke(i, j, k) + m_fuelSmokeRate * deltaFuel;
								m_smoke.set(smoke, i, j, k);
							}
						}

						float temp = m_temperature(i, j, k) + m_fuelTemperatureRate * m_heat(i, j, k);
						m_temperature.set(temp, i, j, k);

						// update divergence as well
						int IX = i * m_pressure.getYSize() * m_pressure.getZSize() + j *  m_pressure.getZSize() + k;
						float divergence = m_divergence(IX) + m_fuelDivergenceRate * deltaFuel / m_timeStep;
						m_divergence(IX) = divergence;
					}
				}
			}
		}
	}
}

void GasSolver::applyFlameCooling()
{
#pragma omp parallel for    
	for (int i = 0; i < m_pressure.getXSize(); i++)
	{
		for (int j = 0; j < m_pressure.getYSize(); j++)
		{
			for (int k = 0; k < m_pressure.getZSize(); k++)
			{
				if (m_layer(i, j, k) != SOLID)
				{
					if (m_heat(i, j, k) > 0)
					{
						float burn = m_heat(i, j, k);
						float controlField = 1.0f;
						if (m_enableFlameCoolingControlField)
							controlField *= smoothstep(m_minBurnThreshold, m_maxBurnThreshold, burn);

						float decay = exp(-m_flameCooling * m_timeStep * controlField) * m_heat(i, j, k);
						m_heat.set(decay, i, j, k);
					}
				}
			}
		}
	}
}

void GasSolver::applyShredding()
{
	float step = 0.5 * m_timeStep * m_fps * m_shreddingAmount;

#pragma omp parallel for    
	for (int i = 0; i < m_pressure.getXSize(); i++)
	{
		for (int j = 0; j < m_pressure.getYSize(); j++)
		{
			for (int k = 0; k < m_pressure.getZSize(); k++)
			{
				if (m_layer(i, j, k) != SOLID)
				{
					vector3 heatgrad = heatGradient(i, j, k);

					if (m_shreddingAmount && heatgrad.length())
					{
						float x = i * m_dx;
						float y = j * m_dx;
						float z = k * m_dx;

						float velu = m_velu(i, j, k) + step * heatgrad.x;
						m_velu.set(velu, i, j, k);

						velu = m_velu(i + 1, j, k) + step * heatgrad.x;
						m_velu.set(velu, i + 1, j, k);

						float velv = m_velv(i, j, k) + step * heatgrad.y;
						m_velv.set(velv, i, j, k);

						velv = m_velv(i, j + 1, k) + step * heatgrad.y;
						m_velv.set(velv, i, j + 1, k);

						float velw = m_velw(i, j, k) + step * heatgrad.z;
						m_velw.set(velw, i, j, k);

						velw = m_velw(i, j, k + 1) + step * heatgrad.z;
						m_velw.set(velw, i, j, k + 1);
					}
				}
			}
		}
	}
}

void GasSolver::autoResizeContainer()
{
}

void GasSolver::setSolveRegion()
{
	int xmax = 0;
	int ymax = 0;
	int zmax = 0;
	int xmin = m_xsize - 1;
	int ymin = m_ysize - 1;
	int zmin = m_zsize - 1;

	// determine bounds of active voxels
#pragma omp parallel for    
	for (int i = 0; i < m_pressure.getXSize(); i++)
	{
		for (int j = 0; j < m_pressure.getYSize(); j++)
		{
			for (int k = 0; k < m_pressure.getZSize(); k++)
			{
				float density = m_smoke(i, j, k);
				float temperature = m_temperature(i, j, k);
				float heat = m_heat(i, j, k);
				float fuel = m_fuel(i, j, k);

				if (density > FIELD_THRESHOLD || temperature > FIELD_THRESHOLD || heat > FIELD_THRESHOLD || fuel > FIELD_THRESHOLD)
				{
					if (i < xmin) xmin = i;
					if (j < ymin) ymin = j;
					if (k < zmin) zmin = k;
					if (i > xmax) xmax = i;
					if (j > ymax) ymax = j;
					if (k > zmax) zmax = k;
				}

				m_layer.set(EMPTY, i, j, k);
			}
		}
	}

	// compute maximum distance density can move in single time step
	float maxv = sqrt(m_velu.maxval() * m_velu.maxval()) + (m_velv.maxval() * m_velv.maxval()) + (m_velw.maxval() * m_velw.maxval());
	float solve_distance = int(m_cflScale * maxv * m_timeStep + 0.5) + m_solveDistance;
	
	// set solve region bounds
	xmin -= solve_distance;
	m_simBounds.xmin = clamp(xmin, 0, m_xsize - 1);
	ymin -= solve_distance;
	m_simBounds.ymin = clamp(ymin, 0, m_ysize - 1);
	zmin -= solve_distance;
	m_simBounds.zmin = clamp(zmin, 0, m_zsize - 1);
	xmax += solve_distance;
	m_simBounds.xmax = clamp(xmax, 0, m_xsize - 1);
	ymax += solve_distance;
	m_simBounds.ymax = clamp(ymax, 0, m_ysize - 1);
	zmax += solve_distance;
	m_simBounds.zmax = clamp(zmax, 0, m_zsize - 1);

#pragma omp parallel for    
	for (int i = m_simBounds.xmin; i <= m_simBounds.xmax; i++)
	{
		for (int j = m_simBounds.ymin; j <= m_simBounds.ymax; j++)
		{
			for (int k = m_simBounds.zmin; k <= m_simBounds.zmax; k++)
			{
				m_layer.set(DENSITY, i, j, k);
			}
		}
	}

	//std::cout << "Solve Distance: " << solve_distance << std::endl;
	//std::cout << "Solve Bounds: " << m_simBounds.xmin << " " << m_simBounds.xmax << " " << m_simBounds.ymin << " " << m_simBounds.ymax << " " << m_simBounds.zmin << " " << m_simBounds.zmax << std::endl;
}