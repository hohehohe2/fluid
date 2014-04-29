#include "CalculatorDensity.h"

#include "Constants.h"
#include <hohe2Common/container/CellCodeCalculator.h>
#include <hohe2Common/container/CompactHash.h>
#include "ParticlesFluid.h"
#include "ParticlesWall.h"
#include "SphKernel.h"


using namespace hohehohe2;


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void CalculatorDensity::calculation_host_(ParticlesFluid& particles, const SphKernel& sphKernel, const CellCodeCalculator& ccc, const CompactHash& cHash,
											   const ParticlesWall* particlesWall, const CompactHash* cHashWall)
{
	particles.m_pos->sync(HOST);
	particles.m_sortedIdMap->sync(HOST);
	const float* pxs = particles.m_pos->xs(HOST);
	const float* pys = particles.m_pos->ys(HOST);
	const float* pzs = particles.m_pos->zs(HOST);
	const unsigned int* sortedIdMaps = particles.m_sortedIdMap->get(HOST);
	float* ds = particles.m_density->get(HOST);

	unsigned int size = particles.size();

	#pragma omp parallel for
	for (int idP = 0; idP < (int)size; ++idP)
	{
		float sumW = 0.0f;

		for (unsigned int i= 0; i < 27; ++i)
		{
			bool isValid;
			unsigned int code = ccc.getNeighborCode32(isValid, pxs[idP], pys[idP], pzs[idP], i);
			if ( ! isValid)
			{
				continue;
			}
			unsigned int index;
			unsigned int numObjects = cHash.lookup(index, code);
			for (unsigned int j = 0; j < numObjects; ++j)
			{
				unsigned int idN = sortedIdMaps[index + j];
				float distx = pxs[idN] - pxs[idP];
				float disty = pys[idN] - pys[idP];
				float distz = pzs[idN] - pzs[idP];
				float dist2 = distx * distx + disty * disty + distz * distz;
				sumW += sphKernel.w(dist2);
			}
		}

		ds[idP] = sumW * m_particleMass;

		//Set rest density as minimum to fake air pressure.
		if (ds[idP] < Constants::RO0)
		{
			ds[idP] = Constants::RO0;
		}

	}

	particles.setClean(HOST);
}
