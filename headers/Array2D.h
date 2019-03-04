//
// Array2D.h 
// Copyright (c) 2016 Pixel Grammar, LLC
// author: Douglas Creel
//

#ifndef ARRAY2D_H
#define ARRAY2D_H

#include <cassert>
#include <iostream>
#include <omp.h>


template <typename T>
class Array2D
{
public:
	Array2D()
	{
	}

	~Array2D()
	{
	}

	void init(int w, int h, T val)
	{
		m_width = w;
		m_height = h;
		grid.resize(w);

#pragma omp parallel for
		for (int i = 0; i < w; i++)
		{
			grid[i].resize(h);
			for (int j = 0; j < h; j++)
			{
				grid[i][j] = val;
			}
		}
	}

	int getWidth()
	{
		return m_width;
	}

	int getHeight()
	{
		return m_height;
	}

	void set(T f, int w, int h)
	{
		grid[w][h] = f;
	}

	typename T& operator() (int w, int h)
	{
		double val = 0;
		if ((w < 0 || w >= m_width) || (h < 0 || h >= m_height)) return val;
		assert(w >= 0 && w < m_width);
		assert(h >= 0 && h < m_height);
		return grid[w][h];
	}

private:
	std::vector< std::vector<T> > grid;
	int m_width, m_height;
};


#endif
