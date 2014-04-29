#ifndef hohe_CalculatorDensity_H
#define hohe_CalculatorDensity_H

#include <hohe2Common/cuda/Buffer.h>

namespace hohehohe2
{

class SphKernel;
struct ParticlesFluid;
class CompactHash;
class CellCodeCalculator;
struct ParticlesWall;

//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Density calculator for SPH.
class CalculatorDensity
{

public:

	///Constructor.
	CalculatorDensity(float particleMass) : m_particleMass(particleMass){}

	///Main method to calculate the acceleration contribution by the pressure force.
	///if particlesWall is given the density is mofied for Akinci2012.
	void calculation(ParticlesFluid& particles, const SphKernel& sphKernel, const CellCodeCalculator& ccc, const CompactHash& cHash, MemoryType mType,
		ParticlesWall* particlesWall=NULL, const CompactHash* cHashWall=NULL)
	{
		if (mType == HOST)
		{
			calculation_host_(particles, sphKernel, ccc, cHash, particlesWall, cHashWall);
		}
		else
		{
			calculation_device_(particles, sphKernel, ccc, cHash, particlesWall, cHashWall);
		}

	}


private:

	///Particle mass.
	float m_particleMass;

private:

	void calculation_host_(ParticlesFluid& particles, const SphKernel& sphKernel, const CellCodeCalculator& ccc, const CompactHash& cHash,
		const ParticlesWall* particlesWall, const CompactHash* cHashWall);
	void calculation_device_(ParticlesFluid& particles, const SphKernel& sphKernel, const CellCodeCalculator& ccc, const CompactHash& cHash,
		const ParticlesWall* particlesWall, const CompactHash* cHashWall)
	{
		//To be implemented.
		assert(false);
	}
};

}

#endif
