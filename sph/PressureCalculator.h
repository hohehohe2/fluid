#ifndef hohe_PressureCalculator_H
#define hohe_PressureCalculator_H

#include <hohe2Common/cuda/Buffer.h>


namespace hohehohe2
{

class SphKernel;
struct FluidParticles;
class CompactHash;


//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Pressure calculator for SPH.
class PressureCalculator
{

public:

	///Constructor.
	PressureCalculator(float particleMass) : m_particleMass(particleMass){}

	///Main method to calculate the acceleration contribution by the pressure force.
	void calcAcceleration(FluidParticles& particles, const SphKernel& sphKernel, MemoryType mType, const CompactHash& cHash)
	{
		if (mType == HOST)
		{
			calcAcceleration_host_(particles, sphKernel, cHash);
		}
		else
		{
			calcAcceleration_device_(particles, sphKernel, cHash);
		}

	}


private:

	///Density->pressure stiffness coefficient.
	/**
	It is not in Constants.h because this is an artificially enough physical value.
	**/
	static const float K;

	///Particle mass.
	float m_particleMass;

private:

	void calcAcceleration_host_(FluidParticles& particles, const SphKernel& sphKernel, const CompactHash& cHash);
	void calcAcceleration_device_(FluidParticles& particles, const SphKernel& sphKernel, const CompactHash& cHash)
	{
		//To be implemented.
		assert(false);
	}

	float densityToPressure_(float density);
};

}

#endif
