//
// PressureSolver.cpp
// Copyright (c) 2016-2018
// author: Douglas Creel
//


#include "PressureSolver.h"



PressureSolver::PressureSolver()
{
}

PressureSolver::PressureSolver(float tol, int maxiter)
{
	m_tolerance = tol;
	m_maxIterations = maxiter;
}

PressureSolver::~PressureSolver()
{
}

void PressureSolver::initMultiGrid(int numGrids, int xsize, int ysize, int zsize)
{
	m_xsize = xsize;
	m_ysize = ysize;
	m_zsize = zsize;

	// compute restriction and interpolation matrices
	m_numGrids = numGrids;
	computeInterpolateMatrices();
	computeRestrictMatrices();
}

VectorXd PressureSolver::conjugateGradientSolver(SparseMatrix<double> coeffs, VectorXd div, VectorXd guess)
{
	// solve Ax = b
	int nvoxels = div.size();
	VectorXd x(nvoxels);

	// Incomplete Cholesky is way slower. What's up?
	//ConjugateGradient<SparseMatrix<double>, Lower|Upper, IncompleteCholesky<double>> solver;
	ConjugateGradient<SparseMatrix<double>, Lower|Upper> solver;

	solver.setMaxIterations(m_maxIterations);
	solver.setTolerance(m_tolerance);

	solver.compute(coeffs);

	if (solver.info() != Success) {
		std::cout << "Decomposition failed!" << std::endl;
	}

	//x = solver.solve(div);
	x = solver.solveWithGuess(div, guess);

	if (solver.info() != Success) {
		std::cout << "Solver did not converge!" << std::endl;
	}

	std::cout << ">>>>>>>> Pressure solve <<<<<<<<" << std::endl;
	std::cout << "#iterations:     " << solver.iterations() << std::endl;
	std::cout << "estimated error: " << solver.error() << std::endl;

	return x;
}

VectorXd PressureSolver::jacobiSolver(bool useGuess, SparseMatrix<double> A, VectorXd b, VectorXd guess)
{
	double maxerror = 0;
	VectorXd xnew, xold;
	int nvoxels = b.size();
	xnew.resize(nvoxels);
	xnew.setZero();
	xold.resize(nvoxels);
	xold.setZero();

	if (useGuess) xold = guess;

	for (int k = 0; k < m_maxIterations; k++)
	{
		for (int l = 0; l < A.outerSize(); ++l)
		{
			double sum = 0;
			for (SparseMatrix<double>::InnerIterator it(A, l); it; ++it)
			{
				if (it.index() != l)
					sum += it.value() * xold(it.index());
			}

			// This is weighted Jacobi
			double weight = 2.0 / 3.0;
			xnew(l) = weight * (b(l) - sum) / A.coeff(l, l) + (1.0 - weight) * xold(l);
	    }

#pragma omp parallel for    
		for (int i = 0; i < nvoxels; i++)
		{
			xold(i) = xnew(i);
		}

		// residual = A * x - b
		VectorXd residual = A * xnew - b;
		maxerror = residual.maxCoeff();
		if (maxerror < m_tolerance)
			break;
    }

	std::cout << "Jacobi solver: max error = " << maxerror << std::endl;

	return xnew;
}

// N-grid V-cycle
VectorXd PressureSolver::multiGridSolver_Vcycle(bool useGuess, SparseMatrix<double> A, VectorXd b, VectorXd guess)
{
	if (!useGuess) guess.setZero();

	int maxiter = m_maxIterations;

	int numSweepsDown = 3;
	int numSweepsUp = 3;

	int nvoxels = m_xsize * m_ysize * m_zsize;

	SparseMatrix<double> An; 
	An.resize(nvoxels, nvoxels);
	An.setZero();
	An = A;
	m_A.push_back(An);

	VectorXd xnew;
	xnew.setZero();

	// relax on fine grid (3 Jacobi)
	m_maxIterations = numSweepsDown;
	xnew = jacobiSolver(true, A, b, guess);

	m_guess.push_back(xnew);

	// compute residual
	VectorXd residual = b - A * xnew;
	//m_rhs.push_back(residual);
	m_rhs.push_back(b);

	VectorXd errorh, errornh;
	guess.setZero();

	std::cout << "DEBUG0 " << std::endl;
	// relax on coarse grids
	for (int i = 0; i < m_numGrids - 1; i++)
	{
		// restrict vectors
		VectorXd residualnh;
		residualnh = restrict(residual, i);
		m_rhs.push_back(residualnh);

		std::cout << "DEBUG1 " << std::endl;
		SparseMatrix<double> A2h;
		A2h.resize(residualnh.size(), residualnh.size());
		std::cout << "DEBUG1A " << m_restrict[i].size() << " " << An.size() << " " << m_interpolate[i].size() << std::endl;
		A2h = m_restrict[i] * An * m_interpolate[i];
		std::cout << "DEBUG2 " << std::endl;

		m_A.push_back(A2h);

		// if coarsest grid then "solve" (should use direct solve instead)
		if (i == m_numGrids - 2)
			m_maxIterations = 100;

		errornh = jacobiSolver(false, A2h, residualnh, guess);
		m_guess.push_back(errornh);

		residual = residualnh - A2h * errornh;

		guess = errornh;

		An.resize(A2h.rows(), A2h.cols());
		An.setZero();
		An = A2h;
	}

	// correct and relax on finer grids
	m_maxIterations = numSweepsUp;
	for (int i = m_numGrids - 2; i >= 0; i--)
	{
		// interpolate to fine grid
		errorh = interpolate(m_guess[i + 1], i);

		// add error back to xnew
		guess = m_guess[i] + errorh;

		// iterate on finer grid (3 Jacobi)
		xnew = jacobiSolver(true, m_A[i], m_rhs[i], guess);
	}

	m_maxIterations = maxiter;

	return xnew;
}

void PressureSolver::setMaxIterations(int iterations)
{
	m_maxIterations = iterations;
}

void PressureSolver::setTolerance(double tolerance)
{
	m_tolerance = tolerance;
}

// 
// private methods
//

void PressureSolver::computeInterpolateMatrices()
{
	int nvoxels = (m_xsize * m_ysize * m_zsize);

	// 27 point kernel
	double F[3][3][3] = {
		{ { 0.125, 0.25, 0.125 },
		{ 0.25, 0.5, 0.25 },
		{ 0.125, 0.25, 0.125 } },
		{ { 0.25, 0.5, 0.25 },
		{ 0.5, 4.0, 0.5 },
		{ 0.25, 0.5, 0.25 } },
		{ { 0.125, 0.25, 0.125 },
		{ 0.25, 0.5, 0.25 },
		{ 0.125, 0.25, 0.125 } }
	};

	for (int idx = 0; idx < m_numGrids; idx++)
	{
		int twohsize = (int)(0.125*nvoxels);

		SparseMatrix<double> interpolate;
		interpolate.resize(nvoxels, twohsize);
		interpolate.setZero();
		interpolate.reserve(VectorXi::Constant(twohsize, 27));

#pragma omp parallel for 
		for (int l = 0; l < twohsize; l++)
		{
			int i = 0; int j = 0; int k = 0;
			int IX = i * m_ysize *  m_zsize + j *  m_zsize + k;
			int IXPI = (i + 1) * m_ysize *  m_zsize + j * m_zsize + k;
			int IXPJ = i * m_ysize *  m_zsize + (j + 1) * m_zsize + k;
			int IXPK = i * m_ysize *  m_zsize + j * m_zsize + (k + 1);
			int IXNI = (i - 1) * m_ysize *  m_zsize + j * m_zsize + k;
			int IXNJ = i * m_ysize *  m_zsize + (j - 1) * m_zsize + k;
			int IXNK = i * m_ysize *  m_zsize + j * m_zsize + (k - 1);

			int IXPIPJ = (i + 1) * m_ysize *  m_zsize + (j + 1) * m_zsize + k;
			int IXPINJ = (i + 1) * m_ysize *  m_zsize + (j - 1) * m_zsize + k;
			int IXPIPK = (i + 1) * m_ysize *  m_zsize + j * m_zsize + (k + 1);
			int IXPINK = (i + 1) * m_ysize *  m_zsize + j * m_zsize + (k - 1);

			int IXNIPJ = (i - 1) * m_ysize *  m_zsize + (j + 1) * m_zsize + k;
			int IXNINJ = (i - 1) * m_ysize *  m_zsize + (j - 1) * m_zsize + k;
			int IXNIPK = (i - 1) * m_ysize *  m_zsize + j * m_zsize + (k + 1);
			int IXNINK = (i - 1) * m_ysize *  m_zsize + j * m_zsize + (k - 1);

			int IXPJPK = i * m_ysize *  m_zsize + (j + 1) * m_zsize + (k + 1);
			int IXNJPK = i * m_ysize *  m_zsize + (j - 1) * m_zsize + (k + 1);
			int IXNJNK = i * m_ysize *  m_zsize + (j - 1) * m_zsize + (k - 1);
			int IXPJNK = i * m_ysize *  m_zsize + (j + 1) * m_zsize + (k - 1);

			int IXPJPI = (i + 1) * m_ysize *  m_zsize + (j + 1) * m_zsize + k;
			int IXPJNI = (i - 1) * m_ysize *  m_zsize + (j + 1) * m_zsize + k;
			int IXNJPI = (i + 1) * m_ysize *  m_zsize + (j - 1) * m_zsize + k;
			int IXNJNI = (i - 1) * m_ysize *  m_zsize + (j - 1) * m_zsize + k;

			int IXPKPI = (i + 1) * m_ysize *  m_zsize + j * m_zsize + (k + 1);
			int IXPKNI = (i - 1) * m_ysize *  m_zsize + j * m_zsize + (k + 1);
			int IXPKPJ = i * m_ysize *  m_zsize + (j + 1) * m_zsize + (k + 1);
			int IXPKNJ = i * m_ysize *  m_zsize + (j - 1) * m_zsize + (k + 1);

			IX += l * 8;
			IXPI += l * 8;
			IXNI += l * 8;
			IXPJ += l * 8;
			IXNJ += l * 8;
			IXPK += l * 8;
			IXNK += l * 8;

			IXPIPJ += l * 8;
			IXPINJ += l * 8;
			IXPIPK += l * 8;
			IXPINK += l * 8;

			IXNIPJ += l * 8;
			IXNINJ += l * 8;
			IXNIPK += l * 8;
			IXNINK += l * 8;

			IXPJPK += l * 8;
			IXNJPK += l * 8;
			IXNJNK += l * 8;
			IXPJNK += l * 8;

			IXPJPI += l * 8;
			IXPJNI += l * 8;
			IXNJPI += l * 8; 
			IXNJNI += l * 8;

			IXPKPI += l * 8;
			IXPKNI += l * 8;
			IXPKPJ += l * 8;
			IXPKNJ += l * 8;

			if (IX < nvoxels && IX >= 0)
				interpolate.coeffRef(IX, l) = F[1][1][1];

			if (IXPI < nvoxels && IXPI >= 0)
				interpolate.coeffRef(IXPI, l) = F[2][1][1];

			if (IXPIPJ < nvoxels && IXPIPJ >= 0)
				interpolate.coeffRef(IXPIPJ, l) = F[2][2][1];
			if (IXPINJ < nvoxels && IXPINJ >= 0)
				interpolate.coeffRef(IXPINJ, l) = F[2][0][1];
			if (IXPIPK < nvoxels && IXPIPK >= 0)
				interpolate.coeffRef(IXPIPK, l) = F[2][1][2];
			if (IXPINK < nvoxels && IXPINK >= 0)
				interpolate.coeffRef(IXPINK, l) = F[2][1][0];

			if (IXNI < nvoxels && IXNI >= 0)
				interpolate.coeffRef(IXNI, l) = F[0][1][1];

			if (IXNIPJ < nvoxels && IXNIPJ >= 0)
				interpolate.coeffRef(IXNIPJ, l) = F[0][2][1];
			if (IXNINJ < nvoxels && IXNINJ >= 0)
				interpolate.coeffRef(IXNINJ, l) = F[0][0][1];
			if (IXNIPK < nvoxels && IXNIPK >= 0)
				interpolate.coeffRef(IXNIPK, l) = F[0][1][2];
			if (IXNINK < nvoxels && IXNINK >= 0)
				interpolate.coeffRef(IXNINK, l) = F[0][1][0];

			if (IXPJ < nvoxels && IXPJ >= 0)
				interpolate.coeffRef(IXPJ, l) = F[1][2][1];

			if (IXPJPI < nvoxels && IXPJPI >= 0)
				interpolate.coeffRef(IXPJPI, l) = F[2][2][1];
			if (IXPJNI < nvoxels && IXPJNI >= 0)
				interpolate.coeffRef(IXPJNI, l) = F[0][2][1];
			if (IXPJPK < nvoxels && IXPJPK >= 0)
				interpolate.coeffRef(IXPJPK, l) = F[1][2][2];
			if (IXPJNK < nvoxels && IXPJNK >= 0)
				interpolate.coeffRef(IXPJNK, l) = F[1][2][0];

			if (IXNJ < nvoxels && IXNJ >= 0)
				interpolate.coeffRef(IXNJ, l) = F[1][0][1];

			if (IXNJPI < nvoxels && IXNJPI >= 0)
				interpolate.coeffRef(IXNJPI, l) = F[2][0][1];
			if (IXNJNI < nvoxels && IXNJNI >= 0)
				interpolate.coeffRef(IXNJNI, l) = F[0][0][1];
			if (IXNJPK < nvoxels && IXNJPK >= 0)
				interpolate.coeffRef(IXNJPK, l) = F[1][0][2];
			if (IXNJNK < nvoxels && IXNJNK >= 0)
				interpolate.coeffRef(IXNJNK, l) = F[1][0][0];

			if (IXPK < nvoxels && IXPK >= 0)
				interpolate.coeffRef(IXPK, l) = F[1][1][2];

			if (IXPKPI < nvoxels && IXPKPI >= 0)
				interpolate.coeffRef(IXPKPI, l) = F[2][1][2];
			if (IXPKNI < nvoxels && IXPKNI >= 0)
				interpolate.coeffRef(IXPKNI, l) = F[0][1][2];
			if (IXPKPJ < nvoxels && IXPKPJ >= 0)
				interpolate.coeffRef(IXPKPJ, l) = F[1][2][2];
			if (IXPKNJ < nvoxels && IXPKNJ >= 0)
				interpolate.coeffRef(IXPKNJ, l) = F[1][0][2];

			if (IXNK < nvoxels && IXNK >= 0)
				interpolate.coeffRef(IXNK, l) = F[1][1][0];
		}

		interpolate.makeCompressed();

		m_interpolate.push_back(interpolate); 

		nvoxels = (int)(0.125 * nvoxels);
	}
}

void PressureSolver::computeRestrictMatrices()
{
#pragma omp parallel for 
	for (int idx = 0; idx < m_numGrids-1; idx++)
	{
		SparseMatrix<double> restrict;

		restrict = 0.125 * m_interpolate[idx].transpose();
		restrict.makeCompressed();
		//std::cout << "DEBUG " << m_interpolate[idx].size() << " " << restrict.size() << std::endl;

		m_restrict.push_back(restrict);
	}
}

VectorXd PressureSolver::restrict(VectorXd residual, int gridNum)
{
	int nvoxels = residual.size();
	int nsize = (int)(0.125 * nvoxels);
	VectorXd residual2h;
	residual2h.resize(nsize);

	// this is just injecting values into new vector from old, should use restriction matrix
//#pragma omp parallel for   
	//for (int i = 0; i < nsize; i++)
	//{
		//residual2h[i] = residual[8 * i];
	//}
	//std::cout << "DEBUG1" << m_restrict[gridNum].rows() << " " << m_restrict[gridNum].cols() << " " << residual.size() << std::endl;
	residual2h = m_restrict[gridNum] * residual;

	return residual2h;
}

VectorXd PressureSolver::interpolate(VectorXd error, int gridNum)
{
	VectorXd error2h;

	error2h = m_interpolate[gridNum] * error;

	return error2h;
}
