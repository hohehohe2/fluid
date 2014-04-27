#include "DensityCalculator.h"

#include "Constants.h"
#include <hohe2Common/container/CellCodeCalculator.h>
#include <hohe2Common/container/CompactHash.h>
#include "FluidParticles.h"
#include "SphKernel.h"


using namespace hohehohe2;


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void DensityCalculator::calcAcceleration_host_(FluidParticles& particles, const SphKernel& sphKernel, const CompactHash& cHash)
{
	particles.m_pos->sync(HOST);
	particles.m_sortedIdMap->sync(HOST);
	const float* pxs = particles.m_pos->xs(HOST);
	const float* pys = particles.m_pos->ys(HOST);
	const float* pzs = particles.m_pos->zs(HOST);
	const unsigned int* sortedIdMaps = particles.m_sortedIdMap->get(HOST);
	float* ds = particles.m_density->get(HOST);

	unsigned int size = particles.size();

	CellCodeCalculator ccc(particles.m_pos->m_lastCalculatedBbox.m_min, sphKernel.r());

	#pragma omp parallel for
	for (int idP = 0; idP < (int)size; ++idP)
	{
		ds[idP] = 0.0f;

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
				ds[idP] += sphKernel.w(dist2);
			}
		}

		ds[idP] *= m_particleMass;

		//Set rest density as minimum to fake air pressure.
		if (ds[idP] < Constants::RO0)
		{
			ds[idP] = Constants::RO0;
		}

	}

	particles.setClean(HOST);
}
