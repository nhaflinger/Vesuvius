//
// PressureSolver.h 
// Copyright (c) 2016-2018
// author: Douglas Creel
//

#ifndef PRESSURESOLVER_H
#define PRESSURESOLVER_H

#include "Vesuvius.h"


class PressureSolver 
{
public:
	PressureSolver();

	PressureSolver(float tol, int maxiter);

	~PressureSolver();

	void initMultiGrid(int numGrids, int xsize, int ysize, int zsize);

	VectorXd conjugateGradientSolver(SparseMatrix<double> coeffs, VectorXd div, VectorXd guess);

	VectorXd jacobiSolver(bool useGuess, SparseMatrix<double> A, VectorXd b, VectorXd guess);

	VectorXd multiGridSolver_Vcycle(bool useGuess, SparseMatrix<double> A, VectorXd b, VectorXd guess);

	void setMaxIterations(int iterations);

	void setTolerance(double tolerance);

private:
	VectorXd restrict(VectorXd residual, int gridNum);

	VectorXd interpolate(VectorXd error, int gridNum);

	void computeInterpolateMatrices();

	void computeRestrictMatrices();

	std::vector<SparseMatrix<double>> m_restrict;

	std::vector<SparseMatrix<double>> m_interpolate;

	std::vector<SparseMatrix<double>> m_A;

	std::vector<VectorXd> m_rhs;

	std::vector<VectorXd> m_guess;

	float m_tolerance = 1.0e-05;

	int m_maxIterations = 100;

	int m_xsize;

	int m_ysize;
	 
	int m_zsize;

	int m_numGrids = 5;
};


#endif
