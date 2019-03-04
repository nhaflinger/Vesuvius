//
// GasSolver.h 
// Copyright (c) 2016-2018
// author: Douglas Creel
//

#ifndef GASSOLVER_H
#define GASSOLVER_H

#include "Vesuvius.h"
#include "Grid3D.h"
#include "PressureSolver.h"
#include "FieldIO.h"
#include "ImplicitSphereSource.h"
#include "ImplicitSphereCollider.h"
#include "GeometrySource.h"
#include "GeometryCollider.h"
#include "PerlinNoise.h"
#include "SimplexNoise.h"
#include "SimplexRDNoise.h"
#include "WaveletNoise.h"
#include "CurlNoise.h"
#include <memory>

using namespace Utilities;


class GasSolver
{
public:
	GasSolver();

	GasSolver(int xsize, int ysize, int zsize, float dx);

	~GasSolver();

	void initialize();

	bool isInitialized();

	void addBodyForce(vector3 force);

	inline int getXSize()
	{
		return m_xsize;
	}

	inline int getYSize()
	{
		return m_ysize;
	}

	inline int getZSize()
	{
		return m_zsize;
	}

	inline void setXSize(int xsize)
	{
		m_xsize = xsize;
	}

	inline void setYSize(int ysize)
	{
		m_ysize = ysize;
	}

	inline void setZSize(int zsize)
	{
		m_zsize = zsize;
	}

	inline void setFrameRate(float fps)
	{
		m_fps = fps;
	}

	inline float getTimeStep()
	{
		return m_timeStep;
	}

	inline float getFrameRate()
	{
		return m_fps;
	}

	inline void setTimeStep(int steps)
	{
		m_subSteps = clamp(steps, 1, 100);
		m_timeStep = (1.0f / m_fps) / m_subSteps;
	}

	inline int getSubSteps()
	{
		return m_subSteps;
	}

	inline void setMaxSubSteps(int numsteps)
	{
		m_maxSubSteps = clamp(numsteps, 1, 100);
	}

	inline void setAdvectionMode(ADVECTMODE mode)
	{
		m_advectMode = mode;
	}

	inline void setTraceMode(TRACEMODE mode)
	{
		m_traceMode = mode;
	}

	void update();

	void saveState(const char* filename);

	void implicitSphereSource(UPDATEBEHAVIOR behavior, bool emptyInterior, float featherAmount, vector3 center, float radius, float density, float densityScale, float temperature, float temperatureScale, float fuel, float fuelScale, bool enableNoise, NoiseParams* noiseparam, SimplexNoise noise, vector3 velocity);

	void implicitSphereCollider(vector3 center, float radius);

	void geometrySource(UPDATEBEHAVIOR behavior, bool emptyInterior, float featherAmount, const char* geometryFile, float density, float densityScale, float temperature, float temperatureScale, float fuel, float fuelScale, bool enableNoise, NoiseParams* noiseparam, SimplexNoise noise, Xform transform, vector3 velocity);

	void geometryCollider(const char* geometryFile, Xform transform);

	float computeCFL();

	inline float simTime()
	{
		return m_simTime;
	}

	inline void setAmbientTemperature(float temperature)
	{
		m_ambientTemperature = temperature;
	}

	inline float getAmbientTemperature()
	{
		return m_ambientTemperature;
	}

	inline void setBouyancyLift(float lift)
	{
		m_bouyancyLift = lift;
	}

	inline void setSmokeWeight(float weight)
	{
		m_smokeWeight = weight;
	}

	inline float getBouyancyList()
	{
		return m_bouyancyLift;
	}

	inline float getSmokeWeight()
	{
		return m_smokeWeight;
	}

	inline float getAmbientPressure()
	{
		return m_ambientPresure;
	}

	inline void setAmbientPressure(float pressure)
	{
		m_ambientPresure = pressure;
	}

	inline void enableViscosity(bool enabled)
	{
		m_enableViscosity = enabled;
	}

	inline void setDissipationAmount(float dissipation)
	{
		m_dissipation = dissipation;
		if (m_dissipation < 0) m_dissipation = 0.f;
	}

	inline void enableDissipation(bool enabled)
	{
		m_enableDissipation = enabled;
	}

	inline void setDiffusionAmount(float diffusion)
	{
		m_diffusion = diffusion;
		if (m_diffusion < 0.f) m_diffusion = 0.f;
	}

	inline void setCoolingAmount(float cooling)
	{
		m_cooling = cooling;
		if (m_cooling < 0.f) m_cooling = 0.f;
	}

	inline void setCFLScale(float scale)
	{
		m_cflScale = scale;
	}

	inline void setBoundaryType(BOUNDARYTYPE type)
	{
		m_boundaryType = type;
	}

	inline void enableVorticityConfinement(bool enabled)
	{
		m_enableVorticityConfinement = enabled;
	}

	inline void enableWindForce(bool enable)
	{
		m_enableWindForce = enable;
	}

	inline void enableTurbulence(bool enabled)
	{
		m_enableTurbulence = enabled;
	}

	inline void setTurbulenceStrength(float strength)
	{
		m_turbulence = strength;
	}

	inline void setTurbulenceThreshold(float threshold)
	{
		m_turbulenceThreshold = threshold;
	}

	inline void enableTurbulenceControlField(bool enabled)
	{
		m_enableTurbulenceControlField = enabled;
	}

	inline void setTurbulenceFlow(float flow)
	{
		m_flow = flow;
	}

	inline void setTurbulenceFrequency(vector3 freq)
	{
		m_curlNoise.setFrequency(freq);
	}

	inline void setTurbulenceSettings(float lacunarity, float gain, int octaves)
	{
		m_curlNoise.setTurbulenceSettings(lacunarity, gain, octaves);
	}

	inline void setPressureSolver(PRESSURESOLVER type)
	{
		m_pressureSolver = type;
	}

	inline void enableFire(bool fire)
	{
		m_enableFire = fire;
	}

	inline void setFuelDissipationAmount(float dissipation)
	{
		m_fuelDissipation = dissipation;
	}

	inline void enableFuelDissipation(bool enabled)
	{
		m_enableFuelDissipation = enabled;
	}

	inline void setBurnRate(float rate)
	{
		m_burnRate = rate;
	}

	inline void setIgnitionTemperature(float temperature)
	{
		m_ignitionTemperature = temperature;
	}

	inline void setFuelSmokeRate(float rate)
	{
		m_fuelSmokeRate = rate;
	}

	inline void setFuelBurnRate(float rate)
	{
		m_fuelBurnRate = rate;
	}

	inline void setFuelTemperatureRate(float rate)
	{
		m_fuelTemperatureRate = rate;
	}

	inline void setFuelDivergenceRate(float rate)
	{
		m_fuelDivergenceRate = rate;
	}

	inline void enableFlameCoolingControlField(bool enabled)
	{
		m_enableFlameCoolingControlField = enabled;
	}

	inline void setFlameCoolingRate(float rate)
	{
		m_flameCooling = rate;
	}

	inline void setKWind(float scale)
	{
		m_kwind = scale;
	}

	inline void setWindForce(vector3 windforce)
	{
		m_windForce = windforce;
	}

	inline void setKinematicViscosity(float viscosity)
	{
		m_kinematicViscosity = viscosity;
	}

	inline void setVorticityAmount(float scale)
	{
		m_vorticityScale = scale;
	}

	inline void setMaxPressureSolveIterations(float iter)
	{
		int m_maxPressureSolveIterations = iter;
		solver.setMaxIterations(iter);
	}

	inline void setPressureSolveTolerance(float tol)
	{
		float m_pressureSolveTolerance = tol;
		solver.setTolerance(tol);
    }

	inline void enableFuelAdvection(bool advect)
	{
		m_advectFuel = advect;
	}

	inline void setFuelAdvectionRate(float rate)
	{
		m_fuelAdvectRate = rate;
	}

	inline void setMinBurnThreshold(float threshold)
	{
		m_minBurnThreshold = threshold;
	}

	inline void setMaxBurnThreshold(float threshold)
	{
		m_maxBurnThreshold = threshold;
	}

	inline void setSmokeFromFire(bool enable)
	{
		m_enableFireSmoke = enable;
	}

	inline void setHeatCutoff(float cutoff)
	{
		m_heatCutoff = cutoff;
	}

	inline void enableDenseSmoke(bool enable)
	{
		m_enableDenseSmoke = enable;
	}

	inline void enableShredding(bool shred)
	{
		m_enableShredding = shred;
	}

	inline void setShreddingAmount(float shred)
	{
		m_shreddingAmount = shred;
	}
	
	inline void setVelocityInterpolationMode(INTERPOLATIONMODE mode)
	{
		m_velocityInterpolationMode = mode;
	}

	inline void setDensityInterpolationMode(INTERPOLATIONMODE mode)
	{
		m_densityInterpolationMode = mode;
	}

	inline void enableResizeContainer(bool enable)
	{
		m_autoResizeContainer = enable;
	}

private:

	void m_initializeSimulation();

	void advanceOneStep();

	void advectVelocity();

	void advectScalars();

	vector3 traceParticle(float x, float y, float z, float dt, float offset);

	vector3 bfecc(float x, float y, float z, float dt, float offset);

	vector3 modifiedMacCormack(float x, float y, float z, float dt, float offset);

	vector3 getVelocity(float x, float y, float z, float offset);

	float getSmoke(float x, float y, float z);

	float getTemperature(float x, float y, float z);

	float getFuel(float x, float y, float z);

	float getHeat(float x, float y, float z);

	vector3 getRest(float x, float y, float z);

	void applyBouyancyForce();

	void applyViscosity();

	void project();

	void computePressureCoefficients();

	void computeNegativeDivergence();

	void pressureUpdate();

	float divergence(int ix, int iy, int iz);

	vector3 pressureGradient(int ix, int iy, int iz);

	vector3 temperatureGradient(int ix, int iy, int iz);

	vector3 heatGradient(int ix, int iy, int iz);

	void applyBoundaryConditions();

	vector3 cellCenterVelocity(int ix, int iy, int iz);

	vector3 faceCenterVelocity(int face, int ix, int iy, int iz);

	void applyVorticityConfinement();

	void applyWindForce();

	void applyDiffusion();

	void applyDissipation();

	void applyCooling();

	void updateFluidSources();

	void updateColliders();

	void applyTurbulence();

	void applyFuelDissipation();

	void burnFuel();

	void applyFlameCooling();

	void applyShredding();

	void autoResizeContainer();

	void setSolveRegion();

	bool m_initialized = false;

	int m_xsize;

	int m_ysize;

	int m_zsize;

	float m_dx;

	float m_timeStep;

	float m_simTime = 0.f;

	ADVECTMODE m_advectMode = TRACE;

	TRACEMODE m_traceMode = RUNGEKUTTA2;

	VectorXd m_divergence;

	SparseMatrix<double> m_coefficients;

	Grid3D<float> m_pressure;

	VectorXd m_guess;

	Grid3D<float> m_velu;

	Grid3D<float> m_velv;

	Grid3D<float> m_velw;

	Grid3D<float> m_smoke;

	Grid3D<float> m_temperature;

	Grid3D<float> m_density;

	Grid3D<LAYER> m_layer;

	Grid3D<float> m_volfrac;

	Grid3D<float> m_solidu;

	Grid3D<float> m_solidv;

	Grid3D<float> m_solidw;

	Grid3D<float> m_restu;

	Grid3D<float> m_restv;

	Grid3D<float> m_restw;

	Grid3D<float> m_fuel;

	Grid3D<float> m_heat;

	PressureSolver solver;

	std::vector<ImplicitSphereSource> m_implicitSphereSources;

	std::vector<GeometrySource> m_geometrySources;

	std::vector<ImplicitSphereCollider> m_implicitSphereColliders;

	std::vector<GeometryCollider> m_geometryColliders;

	vector3 m_gravity = { 0.0f, -9.81f, 0.0f };

	float m_ambientTemperature = 0.0f;

	float m_ambientPresure = 1.0f;

	float m_bouyancyLift = 0.5f;

	float m_smokeWeight = 0.0f;

	float m_fps = 24.0f;

	int m_subSteps = 1;

	bool m_enableVorticityConfinement = false;

	bool m_enableViscosity = false;

	float m_kinematicViscosity = 0.0f;

	float m_vorticityScale = 0.0f;

	float m_diffusion = 0.0f;

	float m_dissipation = 0.0f;

	bool m_enableDissipation = false;

	float m_cooling = 0.0f;

	bool m_enableFire = false;

	float m_ignitionTemperature = 0.0f;

	float m_burnRate = 0.0f;

	bool m_enableFuelDissipation = false;

	float m_fuelDissipation = 0.0f;

	float m_fuelSmokeRate = 0.7f;

	float m_fuelBurnRate = 0.9f;

	float m_fuelTemperatureRate = 0.3f;

	float m_fuelDivergenceRate = 1.0f;

	bool m_advectFuel = false;

	float m_fuelAdvectRate = 0.1f;

	bool m_enableFlameCoolingControlField = false;

	float m_flameCooling = 0.0f;

	float m_minBurnThreshold = 0.0f;

	float m_maxBurnThreshold = 3.0f;

	bool m_enableFireSmoke = true;

	float m_heatCutoff = 0.7f;

	bool m_enableDenseSmoke = true;

	bool m_enableShredding = false;

	float m_shreddingAmount = 0.2f;

	float m_cflScale = 1.5f;

	int m_maxSubSteps = 1;

	BOUNDARYTYPE m_boundaryType = OPEN_BOUNDARY;

	bool m_enableWindForce = false;

	float m_kwind = 0.0f;

	vector3 m_windForce = { 0.f, 0.f, 0.f };

	bool m_enableTurbulence = false;

	float m_turbulence = 0.f;

	float m_turbulenceThreshold = 1.e-06;

	CurlNoise m_curlNoise;

	bool m_enableTurbulenceControlField = false;

	float m_flow = 1.0f;

	PRESSURESOLVER m_pressureSolver = PCG;

	int m_maxPressureSolveIterations = 100;

	float m_pressureSolveTolerance = 1.0e-05;

	INTERPOLATIONMODE m_velocityInterpolationMode = TRILINEAR;

	INTERPOLATIONMODE m_densityInterpolationMode = TRILINEAR;

	bool m_autoResizeContainer = false;

	ContainerBounds m_simBounds;

	int m_solveDistance = 16;

	bool m_computePressureCoefficients = true;
};

#endif
