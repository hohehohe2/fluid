#include "CalculatorPressure.h"

#include "Constants.h"
#include <hohe2Common/container/CellCodeCalculator.h>
#include <hohe2Common/container/CompactHash.h>
#include "ParticlesFluid.h"
#include "particlesWall.h"

using namespace hohehohe2;


//-------------------------------------------------------------------
//-------------------------------------------------------------------
const float CalculatorPressure::K = 20000.0f;


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void CalculatorPressure::calculation_host_(ParticlesFluid& particles, float kernelRadius, const CellCodeCalculator& ccc, const CompactHash& cHash,
										   const ParticlesWall* particlesWall, const CompactHash* cHashWall)
{
	m_sphKernelSpiky.setKernelRadius(kernelRadius);
	m_sphKernelPoly6.setKernelRadius(kernelRadius);

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

	const unsigned int size = particles.size();

	#pragma omp parallel for
	for (int idP = 0; idP < (int)size; ++idP)
	{
		const float densityP = ds[idP];
		const float pressureP = densityToPressure_(densityP);

		//----Pressure Muller 2003.
		//for (unsigned int i= 0; i < 27; ++i)
		//{
		//	bool isFilled;
		//	unsigned int code = ccc.getNeighborCode32(isFilled, pxs[idP], pys[idP], pzs[idP], i);
		//	if ( ! isFilled)
		//	{
		//		continue;
		//	}

		//	unsigned int index;
		//	unsigned int numObjects = m_cHash->lookup(index, code);
		//	for (unsigned int j = 0; j < numObjects; ++j)
		//	{
		//		unsigned int idN = sortedIdMaps[index + j];

		//		Point gradW;
		//		m_sphKernelSpiky.gradWPart(gradW, Point(pxs[idP], pys[idP], pzs[idP]), Point(pxs[idN], pys[idN], pzs[idN]));
		//		gradW *= m_sphKernelSpiky.getConstant(); Do it outside the loop if you uncommnet this code for optimisation.

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
		//{
		//	Point sumGradW(0.0f, 0.0f, 0.0f);

		//	for (unsigned int i= 0; i < 27; ++i)
		//	{
		//		bool isFilled;
		//		const unsigned int code = ccc.getNeighborCode32(isFilled, pxs[idP], pys[idP], pzs[idP], i);
		//		if ( ! isFilled)
		//		{
		//			continue;
		//		}

		//		unsigned int index;
		//		const unsigned int numObjects = cHash.lookup(index, code);
		//		for (unsigned int j = 0; j < numObjects; ++j)
		//		{
		//			const unsigned int idN = sortedIdMaps[index + j];
		//			Point gradWPart;
		//			m_sphKernelSpiky.gradWPart(gradWPart, Point(pxs[idP], pys[idP], pzs[idP]), Point(pxs[idN], pys[idN], pzs[idN]));
		//			sumGradW += gradWPart;
		//		}
		//	}
		//	sumGradW *= m_sphKernelSpiky.getConstant();

		//	const float c = - m_particleMass * pressureP / (densityP * densityP);
		//	axs[idP] += c * sumGradW[0];
		//	ays[idP] += c * sumGradW[1];
		//	azs[idP] += c * sumGradW[2];
		//}

		//----Pressure with incompressible approximation and repulsive force (modified Monaghan2000).
		{
			Point deltaA(0.0f, 0.0f, 0.0f); //Delta acceleration.

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

					Point gradW;
					Point pi(pxs[idP], pys[idP], pzs[idP]);
					Point pj(pxs[idN], pys[idN], pzs[idN]);
					m_sphKernelSpiky.gradWPart(gradW, pi, pj);
					gradW *= m_sphKernelSpiky.getConstant();

					const float initContrib = pressureP / (densityP * densityP); //Coefficiant of pressure acceleration, being canceled by the repulsive force.

					const float DP = 0.8f * kernelRadius;
					const float fAB = m_sphKernelPoly6.wPart((pi - pj).squaredNorm()) / m_sphKernelPoly6.wPart(DP * DP);
					const float N = 0.56f; //Fixed the value with try&error using 2 particles example placement.

					float repulsiveContrib = 0.0f;
					if (pressureP < 0.0f)
					{
						//EPSIRON * initContrib is the amount of cancelation of initContirb when distance to the neighbor is DP, i.e. fAB == 1.
						const float EPSIRON = 0.2f;
						repulsiveContrib = - EPSIRON * initContrib * pow(fAB, N);
					}


					deltaA += (initContrib + repulsiveContrib) * gradW;

				}
			}

			deltaA *= - m_particleMass;
			axs[idP] += deltaA[0];
			ays[idP] += deltaA[1];
			azs[idP] += deltaA[2];
		}

		//----Wall pressure akinci2012.

		Point deltaP(0.0f, 0.0f, 0.0f);
		if (particlesWall)
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
				const unsigned int numObjects = cHashWall->lookup(index, code);
				for (unsigned int j = 0; j < numObjects; ++j)
				{
					const unsigned int idW = wsortedIdMaps[index + j];
					Point gradWPart;
					m_sphKernelSpiky.gradWPart(gradWPart, Point(pxs[idP], pys[idP], pzs[idP]), Point(wpxs[idW], wpys[idW], wpzs[idW]));
					deltaP += wvs[idW] * gradWPart;
				}
			}
		}

		deltaP *= - Constants::RO0 * pressureP / (densityP * densityP) * m_sphKernelSpiky.getConstant();
		axs[idP] += deltaP[0];
		ays[idP] += deltaP[1];
		azs[idP] += deltaP[2];
	}

	particles.setClean(HOST);
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
float CalculatorPressure::densityToPressure_(float density)
{
	return K * (density * density / (Constants::RO0 * Constants::RO0) - 1.0f);
}
