#ifndef hohe_CalculatorIntegrateSemiImplicitEuler_H
#define hohe_CalculatorIntegrateSemiImplicitEuler_H

#include <hohe2Common/cuda/Buffer.h>

namespace hohehohe2
{

class SphKernel;
struct ParticlesSph;
class CompactHash;


//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Semi implicit euler time integration.
class CalculatorIntegrateSemiImplicitEuler
{

public:

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
