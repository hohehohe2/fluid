#include "CalculatorPressure.h"

#include "Constants.h"
#include <hohe2Common/container/CellCodeCalculator.h>
#include <hohe2Common/container/CompactHash.h>
#include "ParticlesFluid.h"
#include "SphKernel.h"


using namespace hohehohe2;


//-------------------------------------------------------------------
//-------------------------------------------------------------------
const float CalculatorPressure::K = 20000.0f;


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void CalculatorPressure::calculation_host_(ParticlesFluid& particles, const SphKernel& sphKernel, const CellCodeCalculator& ccc, const CompactHash& cHash, bool isWall)
{
	particles.m_pos->sync(HOST);
	particles.m_density->sync(HOST);
	particles.m_sortedIdMap->sync(HOST);
	const float* pxs = particles.m_pos->xs(HOST);
	const float* pys = particles.m_pos->ys(HOST);
	const float* pzs = particles.m_pos->zs(HOST);
	const float* ds = particles.m_density->get(HOST);
	const unsigned int* sortedIdMaps = particles.m_sortedIdMap->get(HOST);
	float* axs = particles.m_acceleration->xs(HOST);
	float* ays = particles.m_acceleration->ys(HOST);
	float* azs = particles.m_acceleration->zs(HOST);

	unsigned int size = particles.size();

	#pragma omp parallel for
	for (int idP = 0; idP < (int)size; ++idP)
	{
		float densityP = ds[idP];
		float pressureP = densityToPressure_(densityP);

		//----Pressure Muller 2003.
		//for (unsigned int i= 0; i < 27; ++i)
		//{
		//	bool isValid;
		//	unsigned int code = ccc.getNeighborCode32(isValid, pxs[idP], pys[idP], pzs[idP], i);
		//	if ( ! isValid)
		//	{
		//		continue;
		//	}

		//	unsigned int index;
		//	unsigned int numObjects = m_cHash->lookup(index, code);
		//	for (unsigned int j = 0; j < numObjects; ++j)
		//	{
		//		unsigned int idN = sortedIdMaps[index + j];

		//		float gradW[3];
		//		sphKernel.gradW(gradW, pxs[idP], pys[idP], pzs[idP], pxs[idN], pys[idN], pzs[idN]);

		//		float densityN = ds[idN];
		//		float pressureN = densityToPressure_(densityN);
		//		float c = m_particleMass * (pressureP + pressureN) / (2.0f * densityN);
		//		float gradPx = c * gradW[0];
		//		float gradPy = c * gradW[1];
		//		float gradPz = c * gradW[2];
		//		axs[idP] += - gradPx / densityP;
		//		ays[idP] += - gradPy / densityP;
		//		azs[idP] += - gradPz / densityP;
		//	}
		//}


		//----Pressure with incompressible approximation.
		float sumGradW[3];
		sumGradW[0] = 0.0f;
		sumGradW[1] = 0.0f;
		sumGradW[2] = 0.0f;
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
				float gradW[3];
				sphKernel.gradW(gradW, pxs[idP], pys[idP], pzs[idP], pxs[idN], pys[idN], pzs[idN]);
				sumGradW[0] += gradW[0];
				sumGradW[1] += gradW[1];
				sumGradW[2] += gradW[2];
			}
		}

		float c = - m_particleMass * pressureP / (densityP * densityP);
		axs[idP] += c * sumGradW[0];
		ays[idP] += c * sumGradW[1];
		azs[idP] += c * sumGradW[2];

	}

	particles.setClean(HOST);
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
float CalculatorPressure::densityToPressure_(float density)
{
	return K * (density * density / (Constants::RO0 * Constants::RO0) - 1.0f);
}
