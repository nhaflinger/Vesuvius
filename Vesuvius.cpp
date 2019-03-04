//
// Vesuvius.cpp 
// Copyright (c) 2016-2018
// author: Douglas Creel
//

#include "GasSolver.h"


void
usage(const char* program)
{
	std::cout <<
		"Usage: " << program << " [options]\n" <<
		"Description: Gas solver simulator\n" <<
		"Options:\n" <<
		"    -s            start frame\n" <<
		"    -e            end frame\n" <<
		"    -o            output file name\n" <<
		"    -d            print debugging info\n" << std::endl;
}

int main(int argc, const char * argv[])
{
	int startFrame, endFrame; 

	startFrame = 1;
	endFrame = 50;

	std::cout << "Starting GasSolver" << std::endl;

	const char* program = argv[0];
	std::string outfile;
	bool printDebugInfo = false;

	if (argc == 1)
	{
		usage(program);
		exit(0);
	}

	// parse command line
	for (int n = 1; n < argc; ++n)
	{
		std::string str1(argv[n]);
		if (str1 == "-o")
		{
			if (n < argc - 1)
			{
				std::string str2(argv[n + 1]);
				outfile.assign(str2);
			}
		}
		else if (str1 == "-s")
		{
			startFrame = atoi(argv[n + 1]);
		}
		else if (str1 == "-e")
		{
			endFrame = atoi(argv[n + 1]);
		}
		else if (str1 == "-d")
		{
			printDebugInfo = true;
		}
		else if (str1 == "-h" || str1 == "--help")
		{
			usage(program);
			exit(0);
		}
	}

	// initialize gas solver
	int xsize = 150;
	int ysize = 200;
	int zsize = 150;
	float dx = .1;

	// set number of threads
	omp_set_dynamic(1);
	int numthreads = omp_get_num_procs();
	omp_set_num_threads(numthreads);

	GasSolver gassolver(xsize, ysize, zsize, dx);
	gassolver.initialize();
	std::cout << ">>>>>>>> Initializing domain size (" << gassolver.getXSize() << "," << gassolver.getYSize() << "," << gassolver.getZSize() << ")" << " <<<<<<<<" << std::endl;

	int substeps = 1;
	gassolver.setTimeStep(substeps);
	gassolver.setMaxSubSteps(1);

	vector3 gravity = { 0.0f, -9.81f, 0.0f };
	gassolver.addBodyForce(gravity);

	gassolver.setBouyancyLift(0.5f);

	//gassolver.enableWindForce(true);
	//gassolver.setKWind(0.3f);
	//vector3 winddir = { 2.0f, 0.0f, 0.0f };
	//gassolver.setWindForce(winddir);

	gassolver.enableTurbulence(true);
	gassolver.setTurbulenceStrength(0.2f);
	gassolver.setTurbulenceThreshold(1.e-06);
	gassolver.setTurbulenceFlow(1.0);
	vector3 turbfreq = { 1.0f, 1.0f, 1.0f };
	gassolver.setTurbulenceFrequency(turbfreq);
	gassolver.setTurbulenceSettings(2.0f, 0.5f, 2);

	gassolver.enableVorticityConfinement(true);
    gassolver.setVorticityAmount(2.0);

	gassolver.setDiffusionAmount(0);

	gassolver.enableDissipation(false);
	gassolver.setDissipationAmount(0.1);

	gassolver.setCoolingAmount(0.2);

	gassolver.enableFire(false);
	gassolver.enableFuelDissipation(false);
	gassolver.setBurnRate(0.4);
	gassolver.setFuelDissipationAmount(0.0);
	gassolver.setFuelSmokeRate(0.0);
	gassolver.setFuelTemperatureRate(0.3);
	gassolver.setFuelBurnRate(0.9);
	gassolver.setHeatCutoff(0.3);
	gassolver.setFlameCoolingRate(0.7);

	gassolver.enableShredding(false);
	gassolver.setShreddingAmount(0.3);

	// setup some test conditions
	SimplexNoise noise;
	NoiseParams* noiseParams = new NoiseParams;
	noiseParams->setMultAdd = ADDVALUE;
	noiseParams->noiseAmpScale = 1.0f;
	noiseParams->noiseFreqScale = { 1.0f, 1.0f, 1.0f };
	noiseParams->noiseFlow = 1.0f;
	noiseParams->lacunarity = 2.0f;
	noiseParams->gain = 0.5f;
	noiseParams->octaves = 2;

	vector3 velocity = { 0, 1, 0 };
	vector3 emitcenter = { 7.5f, 2.0f, 7.5f };
    gassolver.implicitSphereSource(ADD, false, 0.1f, emitcenter, 2.0f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f, true, noiseParams, noise, velocity);
	Xform transform; transform.scale = { 10, 10, 10 }; transform.translate = { -75, -20, -75}; transform.rotate = { 0, 0, 0 };
	std::string geometryFile = "C:/Users/mrsunshine/Documents/Visual Studio 2015/Projects/Vesuvius/x64/Release/test/torus_no_blosc.0001.vdb";
	//gassolver.geometrySource(ADD, false, 0.2f, geometryFile.c_str(), 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f, false, noiseParams, noise, transform, velocity);

	//vector3 collidercenter = { 10.0, 8.0, 10.0 };
	//gassolver.implicitSphereCollider(collidercenter, 2.0);

	//transform.translate = { -100, -80, -100 };
	//gassolver.geometryCollider(geometryFile.c_str(), transform);

	gassolver.setAdvectionMode(BFECC);
	gassolver.setTraceMode(RUNGEKUTTA2);
	gassolver.setBoundaryType(OPEN_BOUNDARY);
	//gassolver.setPressureSolver(MULTIGRID);
	gassolver.setVelocityInterpolationMode(TRICUBIC);
	gassolver.setDensityInterpolationMode(TRICUBIC);
	gassolver.enableResizeContainer(false);
	gassolver.setMaxPressureSolveIterations(200);
	gassolver.setPressureSolveTolerance(1.e-05);

	Timer simTimer;
	for (int curframe = startFrame; curframe<=endFrame; curframe++)
	{
		simTimer.setTimer();
		std::cout << ">>>>>>>> Frame = " << curframe << " <<<<<<<<" << std::endl;

		gassolver.update();

		float elapsed = simTimer.elapsedTime();
		std::cout << ">>>>>>>> Time in update() = " << elapsed << " seconds <<<<<<<<" << std::endl;

		char buffer[2048];
		sprintf_s(buffer, outfile.c_str(), curframe);
		gassolver.saveState(buffer);

		// animate noise on source
		noiseParams->noiseFlow = 1.0 * gassolver.simTime();
	}

	return 0;
}





