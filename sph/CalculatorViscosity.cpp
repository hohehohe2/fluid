#include "CalculatorViscosity.h"

#include "Constants.h"
#include <hohe2Common/container/CellCodeCalculator.h>
#include <hohe2Common/container/CompactHash.h>
#include "ParticlesFluid.h"
#include "SphKernel.h"


using namespace hohehohe2;


//-------------------------------------------------------------------
//-------------------------------------------------------------------
const float CalculatorViscosity::MU = 0.1f;


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void CalculatorViscosity::calculation_host_(ParticlesFluid& particles, const SphKernel& sphKernel, const CellCodeCalculator& ccc, const CompactHash& cHash)
{
	particles.m_pos->sync(HOST);
	particles.m_velocity->sync(HOST);
	particles.m_density->sync(HOST);
	particles.m_sortedIdMap->sync(HOST);
	const float* pxs = particles.m_pos->xs(HOST);
	const float* pys = particles.m_pos->ys(HOST);
	const float* pzs = particles.m_pos->zs(HOST);
	const float* vxs = particles.m_velocity->xs(HOST);
	const float* vys = particles.m_velocity->ys(HOST);
	const float* vzs = particles.m_velocity->zs(HOST);
	const float* ds = particles.m_density->get(HOST);
	const unsigned int* sortedIdMaps = particles.m_sortedIdMap->get(HOST);
	float* axs = particles.m_acceleration->xs(HOST);
	float* ays = particles.m_acceleration->ys(HOST);
	float* azs = particles.m_acceleration->zs(HOST);

	const unsigned int size = particles.size();

	#pragma omp parallel for
	for (int idP = 0; idP < (int)size; ++idP)
	{
		for (unsigned int i= 0; i < 27; ++i)
		{
			bool isFilled;
			const unsigned int code = ccc.getNeighborCode32(isFilled, pxs[idP], pys[idP], pzs[idP], i);
			if ( ! isFilled)
			{
				continue;
			}

			unsigned int index;
			const unsigned int numObjects = cHash.lookup(index, code);
			for (unsigned int j = 0; j < numObjects; ++j)
			{
				const unsigned int idN = sortedIdMaps[index + j];

				const float distx = pxs[idN] - pxs[idP];
				const float disty = pys[idN] - pys[idP];
				const float distz = pzs[idN] - pzs[idP];
				const float dist2 = distx * distx + disty * disty + distz * distz;
				const float laplace = sphKernel.laplaceW(dist2);
				const float viscosCoef = MU / ds[idN] * laplace;

                axs[idP] += viscosCoef * (vxs[idN] - vxs[idP]) / m_particleMass;
                ays[idP] += viscosCoef * (vys[idN] - vys[idP]) / m_particleMass;
                azs[idP] += viscosCoef * (vzs[idN] - vzs[idP]) / m_particleMass;
			}
		}

	}

	particles.setClean(HOST);
}
