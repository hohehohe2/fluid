#ifndef hohe_CalculatorVolume_H
#define hohe_CalculatorVolume_H

#include <hohe2Common/cuda/Buffer.h>

namespace hohehohe2
{

class SphKernel;
struct ParticlesWall;
class CompactHash;
class CellCodeCalculator;

//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Volume calculator of SPH wall particles (for Versatile Rigid-Fluid Coupling for Incompressible SPH).
class CalculatorVolume
{

public:

	///Main method to calculate the acceleration contribution by the pressure force.
	void calculation(ParticlesWall& particles, const SphKernel& sphKernel, const CellCodeCalculator& ccc, const CompactHash& cHash, MemoryType mType)
	{
		if (mType == HOST)
		{
			calculation_host_(particles, sphKernel, ccc, cHash);
		}
		else
		{
			calculation_device_(particles, sphKernel, ccc, cHash);
		}

	}

private:

	void calculation_host_(ParticlesWall& particles, const SphKernel& sphKernel, const CellCodeCalculator& ccc, const CompactHash& cHash);
	void calculation_device_(ParticlesWall& particles, const SphKernel& sphKernel, const CellCodeCalculator& ccc, const CompactHash& cHash)
	{
		//To be implemented.
		assert(false);
	}
};

}

#endif
