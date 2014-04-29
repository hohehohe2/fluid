#ifndef hohe_CalculatorPressure_H
#define hohe_CalculatorPressure_H

#include <hohe2Common/cuda/Buffer.h>


namespace hohehohe2
{

class SphKernel;
struct ParticlesFluid;
struct ParticlesWall;
class CompactHash;
class CellCodeCalculator;

//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Pressure calculator for SPH.
class CalculatorPressure
{

public:

	///Constructor.
	CalculatorPressure(float particleMass) : m_particleMass(particleMass){}

	///Main method to calculate the acceleration contribution by the pressure force.
	void calculation(ParticlesFluid& particles, const SphKernel& sphKernel, const CellCodeCalculator& ccc, const CompactHash& cHash, MemoryType mType,
		const ParticlesWall* particlesWall=NULL, const CompactHash* cHashWall=NULL)
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

	///Density->pressure stiffness coefficient.
	/**
	It is not in Constants.h because this is an artificially enough physical value.
	**/
	static const float K;

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

	float densityToPressure_(float density);
};

}

#endif
