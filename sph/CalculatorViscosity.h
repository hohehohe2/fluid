#ifndef hohe_CalculatorViscosity_H
#define hohe_CalculatorViscosity_H

#include <hohe2Common/cuda/Buffer.h>

namespace hohehohe2
{

class SphKernel;
struct ParticlesFluid;
class CompactHash;
class CellCodeCalculator;

//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Pressure calculator for SPH.
class CalculatorViscosity
{

public:

	///Constructor.
	CalculatorViscosity(float particleMass) : m_particleMass(particleMass){}

	///Main method to calculate the acceleration contribution by the pressure force.
	void calcAcceleration(ParticlesFluid& particles, const SphKernel& sphKernel, const CellCodeCalculator& ccc, const CompactHash& cHash, MemoryType mType)
	{
		if (mType == HOST)
		{
			calcAcceleration_host_(particles, sphKernel, ccc, cHash);
		}
		else
		{
			calcAcceleration_device_(particles, sphKernel, ccc, cHash);
		}

	}


private:

	///Viscosity coefficient.
	static const float MU;

	///Particle mass.
	float m_particleMass;

private:

	void calcAcceleration_host_(ParticlesFluid& particles, const SphKernel& sphKernel, const CellCodeCalculator& ccc, const CompactHash& cHash);
	void calcAcceleration_device_(ParticlesFluid& particles, const SphKernel& sphKernel, const CellCodeCalculator& ccc, const CompactHash& cHash)
	{
		//To be implemented.
		assert(false);
	}
};

}

#endif
