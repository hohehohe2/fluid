#ifndef hohe_CalculatorViscosity_H
#define hohe_CalculatorViscosity_H

#include <hohe2Common/cuda/Buffer.h>
#include "kernels/SphKernelViscosity.h"

namespace hohehohe2
{

struct ParticlesFluid;
class CompactHash;
class CellCodeCalculator;
struct GlobalFluidParameters;

//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Viscosity force (acceleration) calculator for SPH.
class CalculatorViscosity
{

public:

	///Constructor.
	CalculatorViscosity(float particleMass) : m_particleMass(particleMass){}

	///Main method to calculate the acceleration contribution by the pressure force.
	void calculation(ParticlesFluid& particles, const GlobalFluidParameters& globalParam, const CellCodeCalculator& ccc, const CompactHash& cHash, MemoryType mType)
	{
		if (mType == HOST)
		{
			calculation_host_(particles, globalParam, ccc, cHash);
		}
		else
		{
			calculation_device_(particles, globalParam, ccc, cHash);
		}

	}


private:

	///Kernel for viscosity calculation.
	SphKernelViscosity m_sphKernelViscosity;

	///Viscosity coefficient.
	static const float MU;

	///Particle mass.
	const float m_particleMass;

private:

	void calculation_host_(ParticlesFluid& particles, const GlobalFluidParameters& globalParam, const CellCodeCalculator& ccc, const CompactHash& cHash);
	void calculation_device_(ParticlesFluid& particles, const GlobalFluidParameters& globalParam, const CellCodeCalculator& ccc, const CompactHash& cHash)
	{
		//To be implemented.
		assert(false);
	}
};

}

#endif
