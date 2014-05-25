#include "CalculatorDensity.h"

#include "Constants.h"
#include <hohe2Common/container/CellCodeCalculator.h>
#include <hohe2Common/container/CompactHash.h>
#include "ParticlesFluid.h"
#include "ParticlesWall.h"
#include "FluidSolverSimpleSph.h"


using namespace hohehohe2;


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void CalculatorDensity::calculation_host_(ParticlesFluid& particles, const GlobalFluidParameters& globalParam, const CellCodeCalculator& ccc, const CompactHash& cHash,
										  const ParticlesWall* particlesWall, const CompactHash* cHashWall)
{
	m_sphKernelPoly6.setKernelRadius(globalParam.m_kernelRadius);

	particles.m_pos->sync(HOST);
	particles.m_sortedIdMap->sync(HOST);
	const float* pxs = particles.m_pos->xs(HOST);
	const float* pys = particles.m_pos->ys(HOST);
	const float* pzs = particles.m_pos->zs(HOST);
	const unsigned int* sortedIdMaps = particles.m_sortedIdMap->get(HOST);
	float* ds = particles.m_density->get(HOST);

	const float* wpxs;
	const float* wpys;
	const float* wpzs;
	const float* wvs;
	const unsigned int* wsortedIdMaps;
	if (particlesWall)
	{
		particlesWall->m_pos->sync(HOST);
		particlesWall->m_volume->sync(HOST);
		particlesWall->m_sortedIdMap->sync(HOST);
		wpxs = particlesWall->m_pos->xs(HOST);
		wpys = particlesWall->m_pos->ys(HOST);
		wpzs = particlesWall->m_pos->zs(HOST);
		wvs = particlesWall->m_volume->get(HOST);
		wsortedIdMaps = particlesWall->m_sortedIdMap->get(HOST);
	}

	unsigned int size = particles.size();

	#pragma omp parallel for
	for (int idP = 0; idP < (int)size; ++idP)
	{
		float sumW = 0.0f;

		//----Contribution from neightbor fruid particles.
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
				sumW += m_sphKernelPoly6.wPart(dist2);
			}
		}
		sumW *= m_sphKernelPoly6.getConstant();

		ds[idP] = sumW * m_particleMass;


		//----Contribution from neightbor wall particles [Akinci2012].
		if (particlesWall)
		{
			float densContribWall = 0.0f;
			for (unsigned int i= 0; i < 27; ++i)
			{
				bool isFilled;
				const unsigned int code = ccc.getNeighborCode32(isFilled, pxs[idP], pys[idP], pzs[idP], i);
				if ( ! isFilled)
				{
					continue;
				}
				unsigned int index;
				const unsigned int numObjects = cHashWall->lookup(index, code);
				for (unsigned int j = 0; j < numObjects; ++j)
				{
					const unsigned int idW = wsortedIdMaps[index + j];
					const float distx = wpxs[idW] - pxs[idP];
					const float disty = wpys[idW] - pys[idP];
					const float distz = wpzs[idW] - pzs[idP];
					const float dist2 = distx * distx + disty * disty + distz * distz;
					densContribWall += Constants::RO0 * wvs[idW] * m_sphKernelPoly6.wPart(dist2);
				}
			}
			densContribWall *= m_sphKernelPoly6.getConstant();

			ds[idP] += densContribWall;
		}

		////----Set rest density as minimum to fake air pressure.
		//if (ds[idP] < Constants::RO0)
		//{
		//	ds[idP] = Constants::RO0;
		//}

	}

	particles.setClean(HOST);
}
