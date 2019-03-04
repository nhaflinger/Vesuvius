//
// Utilities.h 
// Copyright (c) 2016-2018
// author: Douglas Creel
//

#ifndef UTILITIES_H
#define UTILITIES_H


#include <math.h>
#include <algorithm> 


#ifdef WIN32
#undef min
#undef max
#endif


using std::min;
using std::max;


struct vector3 {
	float x;
	float y;
	float z;

	vector3 operator=(vector3 a) {
		x = a.x;
		y = a.y;
		z = a.z;
		return *this;
	}

	vector3 operator*(float a) {
		x *= a;
		y *= a;
		z *= a;
		return *this;
	}

	vector3 operator*(vector3 a) {
		x *= a.x;
		y *= a.y;
		z *= a.z;
		return *this;
	}

	vector3 operator+(vector3 a) {
		x += a.x;
		y += a.y;
		z += a.z;
		return *this;
	}

	float length()
	{
		return sqrt(x*x + y*y + z*z);
	}
};


namespace Utilities
{
	inline extern float clamp(float x, float low, float high)
	{
		if (x > high) return high;
		if (x < low) return low;
		else return x;
	}

	inline extern float smoothstep(float edge0, float edge1, float x)
	{
		x = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);

		return x*x*(3 - 2 * x);
	}

    inline extern float smootherstep(float edge0, float edge1, float x)
	{
		x = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);

		return x*x*x*(x*(x * 6 - 15) + 10);
	}

	inline float lerp(float v0, float v1, float t)
	{
		return (1 - t)*v0 + t*v1;
	}
};

#endif
