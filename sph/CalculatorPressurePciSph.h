#ifndef hohe_CalculatorPressurePciSph_H
#define hohe_CalculatorPressurePciSph_H

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
///PCISPH pressure force (acceleration) calculator.
class CalculatorPressurePciSph
{

public:

	///Constructor. PCISPH precomputation.
    CalculatorPressurePciSph::CalculatorPressurePciSph(float particleMass, const SphKernel& sphKernel, float maxRelativeDensityError, unsigned int numMaxIterations)
		: m_particleMass(particleMass), m_sphKernel(sphKernel), m_deltaT(FLT_MAX), m_maxRelativeDensityError(maxRelativeDensityError), m_numMaxIterations(numMaxIterations),
		  m_lastPrecomputeEquilibriumDistance(FLT_MAX), m_lastPrecomputeKernelRadiusPerEquilibriumDistance(UINT_MAX){}

	///PCISPH precomputation, compute m_delta.
	void precompute(float equilibriumDistance, int kernelRadiusPerEquilibriumDistance, float deltaT);

	///Main method to calculate the acceleration contribution by the pressure force.
	void calculation(ParticlesFluid& particles, const CellCodeCalculator& ccc, const CompactHash& cHash, MemoryType mType)
	{
		if (mType == HOST)
		{
			calculation_host_(particles, ccc, cHash);
		}
		else
		{
			calculation_device_(particles, ccc, cHash);
		}

	}


private:

	///Particle mass.
	const float m_particleMass;

	///SPH kernel to use.
	const SphKernel& m_sphKernel;

	///Time step.
	float m_deltaT;

	///PCISPH scaling factor constant. (density error) * m_delta = (presure change needed to reverse the error).
	float m_delta;

    ///Max allowed density error. One more acceleration update occurs after the density meets this condition.
    float m_maxRelativeDensityError;

    ///Max number of prediction-correction iterations.
    unsigned int m_numMaxIterations;

	///equilibriumDistance parameter value when precompute() was called.
	float m_lastPrecomputeEquilibriumDistance;

	///kernelRadiusPerEquilibriumDistance parameter value when precompute() was called.
	unsigned int m_lastPrecomputeKernelRadiusPerEquilibriumDistance;

private:

	void calculation_host_(ParticlesFluid& particles, const CellCodeCalculator& ccc, const CompactHash& cHash);
	void calculation_device_(ParticlesFluid& particles, const CellCodeCalculator& ccc, const CompactHash& cHash)
	{
		//To be implemented.
		assert(false);
	}

	float densityToPressure_(float density);
};

}

#endif
