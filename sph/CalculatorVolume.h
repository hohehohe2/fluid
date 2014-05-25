#ifndef hohe_CalculatorVolume_H
#define hohe_CalculatorVolume_H

#include <hohe2Common/cuda/Buffer.h>
#include "kernels/SphKernelPoly6.h"

namespace hohehohe2
{

struct ParticlesWall;
class CompactHash;
class CellCodeCalculator;
struct GlobalFluidParameters;

//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Volume calculator of SPH wall particles (for Versatile Rigid-Fluid Coupling for Incompressible SPH).
class CalculatorVolume
{

public:

	///Main method to calculate the acceleration contribution by the pressure force.
	void calculation(ParticlesWall& particles, const GlobalFluidParameters& globalParam, const CellCodeCalculator& ccc, const CompactHash& cHash, MemoryType mType)
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

	SphKernelPoly6 m_sphKernelPoly6;

private:

	void calculation_host_(ParticlesWall& particles, const GlobalFluidParameters& globalParam, const CellCodeCalculator& ccc, const CompactHash& cHash);
	void calculation_device_(ParticlesWall& particles, const GlobalFluidParameters& globalParam, const CellCodeCalculator& ccc, const CompactHash& cHash)
	{
		//To be implemented.
		assert(false);
	}
};

}

#endif
