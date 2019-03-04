//
// Timer.h 
// Copyright (c) 2016-2018
// author: Douglas Creel
//

#ifndef TIMER_H
#define TIMER_H

#include "Vesuvius.h"

class Timer {
public:
	Timer() 
	{ 
		m_start = clock();
	}

	virtual ~Timer() { }

	clock_t setTimer()
	{
		clock_t now = clock();
		m_start = now;
		return now;
	}

	float elapsedTime()
	{
		m_end = clock();
		m_elapsed = m_end - m_start;
		m_start = m_end;
		float elapsed = ((float)m_elapsed) / CLOCKS_PER_SEC;
		return elapsed;
	}

private:
	clock_t m_start;
	clock_t m_end;
	clock_t m_elapsed;
};

#endif
