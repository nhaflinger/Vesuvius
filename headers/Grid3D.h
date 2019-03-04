//
// Grid3D.h 
// Copyright (c) 2016-2018
// author: Douglas Creel
//

#ifndef GRID3D_H
#define GRID3D_H

#include "Vesuvius.h"


enum GRIDBEHAVIOR {
	CLAMPED,
	REPEAT,
	ZERO
};


template <typename T>
class Grid3D
{
public:
	Grid3D()
	{
	}

	~Grid3D()
	{
	}

	void init(int nx, int ny, int nz, T val)
	{
		m_nvoxels = nx*ny*nz;
		m_xsize = nx;
		m_ysize = ny;
		m_zsize = nz;

		data.resize(m_nvoxels); 
#pragma omp parallel for   
		for (int i = 0; i < m_nvoxels; i++)
		{
			data[i] = val;
		}

	}

	int IX(int x, int y, int z)
	{
        if (!(x < m_xsize)) x = m_xsize - 1;
		if (!(y < m_ysize)) y = m_ysize - 1;
		if (!(z < m_zsize)) z = m_zsize - 1;
		if (!(x > 0)) x = 0;
		if (!(y > 0)) y = 0;
		if (!(z > 0)) z = 0;

		return x * m_ysize * m_zsize + y * m_zsize + z;
	}

	inline void set(T f, int x, int y, int z)
	{
		data[IX(x, y, z)] = f;
	}

	typename T& operator() (int x, int y, int z)
	{
		if (m_behavior == CLAMPED)
		{
			if (!(x < m_xsize)) x = m_xsize - 1;
			if (!(y < m_ysize)) y = m_ysize - 1;
			if (!(z < m_zsize)) z = m_zsize - 1;
			if (!(x > 0)) x = 0;
			if (!(y > 0)) y = 0;
			if (!(z > 0)) z = 0;
		}
		else if (m_behavior == ZERO)
		{
			T retval;
			if ((x < 0 || x > m_xsize) || (y < 0 || y > m_ysize) || (z < 0 || z > m_zsize)) return retval;
		}

		assert(x >= 0 && x < m_xsize);
		assert(y >= 0 && y < m_ysize);
		assert(z >= 0 && z < m_zsize);

		return data[IX(x, y, z)];
	}

	inline typename T& operator[] (int i)
	{
		assert(i >= 0 && i < m_nvoxels);
		return data[i];
	}

	inline void setIndexValue(int i, T f)
	{
		data[i] = f;
	}

	inline int getNumVoxels()
	{
		return m_nvoxels;
	}

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

	inline void setBoundaryBehavior(GRIDBEHAVIOR mode)
	{
		m_behavior = mode;
	}

	inline void setZero()
	{
#pragma omp parallel for    
		for (int i = 0; i < m_nvoxels; i++)
		{
			data[i] = 0.0;
		}
	}

	inline T maxval() const
	{
		T r = 0;
		for (int i = 0; i < m_nvoxels; i++)
			if (!(std::fabs(data[i]) <= r)) r = std::fabs(data[i]);
		return r;
	}

	inline void clear()
	{
		data.clear();
	}

	typename T trilinearInterpolation(T x, T y, T z)
	{
		T retval;
		T c0, c1, c2, c3, c4, c5, c6, c7;

		int i = (int)abs(floor(x));
		int j = (int)abs(floor(y));
		int k = (int)abs(floor(z));

		c0 = (i + 1 - x) * (j + 1 - y) * (k + 1 - z) * data[IX(i, j, k)];
		c1 = (x - i) * (j + 1 - y) * (k + 1 - z) * data[IX(i + 1, j, k)];
		c2 = (i + 1 - x) * (y - j) * (k + 1 - z) * data[IX(i, j + 1, k)];
		c3 = (x - i) * (y - j) * (k + 1 - z) * data[IX(i + 1, j + 1, k)];
		c4 = (i + 1 - x) * (j + 1 - y) * (z - k) *  data[IX(i, j, k + 1)];
		c5 = (x - i) * (j + 1 - y) * (z - k) * data[IX(i + 1, j, k + 1)];
		c6 = (i + 1 - x) * (y - j) * (z - k) * data[IX(i, j + 1, k + 1)];
		c7 = (x - i) * (y - j) * (z - k) * data[IX(i + 1, j + 1, k + 1)];

		retval = c0 + c1 + c2 + c3 + c4 + c5 + c6 + c7;

		return retval;
	}

	//
	// code attribution: "Visual Simulation of Smoke", Fedkiw, Stam, Jensen
	//
	typename T monotonicCubicInterpolation(T t, T f[4])
	{
		T retval;

		T d_k = T(0.5) * (f[2] - f[0]);
		T d_k1 = T(0.5) * (f[3] - f[1]);
		T delta_k = f[2] - f[1];

		if (delta_k == static_cast<T>(0))
		{
			d_k = static_cast<T>(0);
			d_k1 = static_cast<T>(0);
		}

		T a0 = f[1];
		T a1 = d_k;
		T a2 = (T(3) * delta_k) - (T(2) * d_k) - d_k1;
		T a3 = d_k + d_k1 - (T(2) * delta_k);

		T t1 = t;
		T t2 = t * t;
		T t3 = t2 * t1;

		retval = a3 * t3 + a2 * t2 + a1 * t1 + a0;

		return retval;
	}

	typename T cubicInterpolation(T x, T p[4])
	{
		T retval;

		retval = p[1] + 0.5 * x * (p[2] - p[0] + x * (2.0 * p[0] - 5.0 * p[1] + 4.0 * p[2] - p[3] + x * (3.0 * (p[1] - p[2]) + p[3] - p[0])));

		return retval;
	}

	typename T bicubicInterpolation(T x, T y, T p[4][4])
	{
		T retval;

		T interps[4];
		interps[0] = cubicInterpolation(y, p[0]);
		interps[1] = cubicInterpolation(y, p[1]);
		interps[2] = cubicInterpolation(y, p[2]);
		interps[3] = cubicInterpolation(y, p[3]);

		//retval = monotonicCubicInterpolation(x, interps);
		retval = cubicInterpolation(x, interps);

		return retval;
	}

	typename T tricubicInterpolation(T x, T y, T z)
	{
		T retval;

		int i = (int)abs(floor(x));
		int j = (int)abs(floor(y));
		int k = (int)abs(floor(z));

		T p[4][4][4]; 

		for (int nj = j - 1, xj = 0; nj < j + 3; nj++, xj++)
		{
			for (int nk = k - 1, xk = 0; nk < k + 3; nk++, xk++)
			{
				p[0][xj][xk] = data[IX(i - 1, nj, nk)];
				p[1][xj][xk] = data[IX(i,     nj, nk)]; 
				p[2][xj][xk] = data[IX(i + 1, nj, nk)];
				p[3][xj][xk] = data[IX(i + 2, nj, nk)];
			}
		}
		
		T interps[4]; 
		interps[0] = bicubicInterpolation(y - j, z - k, p[0]);
		interps[1] = bicubicInterpolation(y - j, z - k, p[1]);
		interps[2] = bicubicInterpolation(y - j, z - k, p[2]);
		interps[3] = bicubicInterpolation(y - j, z - k, p[3]);

		//retval = monotonicCubicInterpolation(x - i, interps);
		retval = cubicInterpolation(x - i, interps);

		return retval;
	}

private:
	std::vector<T> data;
	int m_xsize, m_ysize, m_zsize;
	int m_nvoxels;
	GRIDBEHAVIOR m_behavior = CLAMPED;
};


#endif

