#ifndef hohe_CalculatorSemiImplicitEulerIntegrate_H
#define hohe_CalculatorSemiImplicitEulerIntegrate_H

#include <hohe2Common/cuda/Buffer.h>

namespace hohehohe2
{

class SphKernel;
struct ParticlesFluid;
class CompactHash;


//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Semi implicit euler time integration.
class CalculatorSemiImplicitEulerIntegrate
{

public:

	///Main method to integrate.
	void integrate(ParticlesFluid& particles, float deltaT, MemoryType mType)
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

	void integrate_host_(ParticlesFluid& particles, float deltaT);
	void integrate_device_(ParticlesFluid& particles, float deltaT)
	{
		//To be implemented.
		assert(false);
	}
};

}

#endif
