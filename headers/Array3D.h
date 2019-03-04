//
// Array3D.h 
// Copyright (c) 2016 Pixel Grammar, LLC
// author: Douglas Creel
//

#ifndef ARRAY3D_H
#define ARRAY3D_H

#include <cassert>
#include <iostream>


template <typename T>
class Array3D
{
public:
	Array3D()
	{
	}

	~Array3D()
	{
	}

	void init(int w, int h, int d, T val)
	{
		m_width = w;
		m_height = h;
		m_depth = d;
		grid.resize(w);

#pragma omp parallel for
		for (int i = 0; i < w; i++)
		{
			grid[i].resize(h);
			for (int j = 0; j < h; j++)
			{
				grid[i][j].resize(d);
				for (int k = 0; k < d; k++)
				{
					grid[i][j][k] = val;
				}
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

	int getDepth()
	{
		return m_depth;
	}

	void set(T f, int w, int h, int d)
	{
		grid[w][h][d] = f;
	}

	typename T& operator() (int w, int h, int d)
	{
		double val = 0;
		if ((w < 0 || w >= m_width) || (h < 0 || h >= m_height) || (d < 0 || d >= m_depth)) return val;
		assert(w >= 0 && w < m_width);
		assert(h >= 0 && h < m_height);
		assert(d >= 0 && d < m_depth);
		return grid[w][h][d];
	}

private:
	std::vector<std::vector<std::vector<T> > > grid;
	int m_width, m_height, m_depth;
};


#endif