#ifndef hohe_ViscosityCalculator_H
#define hohe_ViscosityCalculator_H

#include <hohe2Common/cuda/Buffer.h>

namespace hohehohe2
{

class SphKernel;
struct FluidParticles;
class CompactHash;


//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Pressure calculator for SPH.
class ViscosityCalculator
{

public:

	///Constructor.
	ViscosityCalculator(float particleMass) : m_particleMass(particleMass){}

	///Main method to calculate the acceleration contribution by the pressure force.
	void calcAcceleration(FluidParticles& particles, const SphKernel& sphKernel, const CompactHash& cHash, MemoryType mType)
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

	///Viscosity coefficient.
	static const float MU;

	///Particle mass.
	float m_particleMass;

private:

	void calcAcceleration_host_(FluidParticles& particles, const SphKernel& sphKernel, const CompactHash& cHash);
	void calcAcceleration_device_(FluidParticles& particles, const SphKernel& sphKernel, const CompactHash& cHash)
	{
		//To be implemented.
		assert(false);
	}
};

}

#endif
