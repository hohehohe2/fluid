#ifndef hohe_SemiImplicitEulerIntegrateCalculator_H
#define hohe_SemiImplicitEulerIntegrateCalculator_H

#include <hohe2Common/cuda/Buffer.h>

namespace hohehohe2
{

class SphKernel;
struct FluidParticles;
class CompactHash;


//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Semi implicit euler time integration.
class SemiImplicitEulerIntegrateCalculator
{

public:

	///Main method to integrate.
	void integrate(FluidParticles& particles, float deltaT, MemoryType mType)
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

	void integrate_host_(FluidParticles& particles, float deltaT);
	void integrate_device_(FluidParticles& particles, float deltaT)
	{
		//To be implemented.
		assert(false);
	}
};

}

#endif
