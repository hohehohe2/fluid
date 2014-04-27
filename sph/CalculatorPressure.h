#ifndef hohe_CalculatorPressure_H
#define hohe_CalculatorPressure_H

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
class CalculatorPressure
{

public:

	///Constructor.
	CalculatorPressure(float particleMass) : m_particleMass(particleMass){}

	///Main method to calculate the acceleration contribution by the pressure force.
	void calcAcceleration(ParticlesFluid& particles, const SphKernel& sphKernel, const CellCodeCalculator& ccc, const CompactHash& cHash, bool isWall, MemoryType mType)
	{
		if (mType == HOST)
		{
			calcAcceleration_host_(particles, sphKernel, ccc, cHash, isWall);
		}
		else
		{
			calcAcceleration_device_(particles, sphKernel, ccc, cHash, isWall);
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

	void calcAcceleration_host_(ParticlesFluid& particles, const SphKernel& sphKernel, const CellCodeCalculator& ccc, const CompactHash& cHash, bool isWall);
	void calcAcceleration_device_(ParticlesFluid& particles, const SphKernel& sphKernel, const CellCodeCalculator& ccc, const CompactHash& cHash, bool isWall)
	{
		//To be implemented.
		assert(false);
	}

	float densityToPressure_(float density);
};

}

#endif
