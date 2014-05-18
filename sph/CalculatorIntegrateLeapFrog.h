#ifndef hohe_CalculatorIntegrateLeapFrog_H
#define hohe_CalculatorIntegrateLeapFrog_H

#include <hohe2Common/cuda/Buffer.h>
#include "ParticlesSph.h"

namespace hohehohe2
{

class SphKernel;
class CompactHash;


//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Semi implicit euler time integration.
class CalculatorIntegrateLeapFrog
{

public:

	///Copy particle's m_acceleration to m_prevAcceleration. Call this at the beginning of the step.
	void copyAcceleration(ParticlesSph& particles, MemoryType mType)
	{
		particles.m_prevAcceleration->copyAll(*particles.m_acceleration, mType);
	}

	///Main method to integrate.
	void integrate(ParticlesSph& particles, float deltaT, MemoryType mType)
	{
		if (mType == HOST)
		{
			integrate_host_(particles, deltaT);
		}
		else
		{
			integrate_device_(particles, deltaT);
		}
	}

private:

	void integrate_host_(ParticlesSph& particles, float deltaT);
	void integrate_device_(ParticlesSph& particles, float deltaT)
	{
		//To be implemented.
		assert(false);
	}
};

}

#endif
