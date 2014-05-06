#include "CalculatorPressurePciSph.h"

#include "Constants.h"
#include <hohe2Common/container/CellCodeCalculator.h>
#include <hohe2Common/container/CompactHash.h>
#include "ParticlesFluid.h"
#include "particlesWall.h"
#include "SphKernel.h"


using namespace hohehohe2;


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void CalculatorPressurePciSph::precompute(float equilibriumDistance, int kernelRadiusPerEquilibriumDistance)
{
	const unsigned int numLines = kernelRadiusPerEquilibriumDistance * 2 + 1;
	const unsigned int numParticles = numLines * numLines * numLines;
	const unsigned int centerPerticleId = numLines * numLines * kernelRadiusPerEquilibriumDistance + numLines * kernelRadiusPerEquilibriumDistance + kernelRadiusPerEquilibriumDistance;

	//Create prototype particles with uniformly filled placement.
	ParticlesFluid* particles = ParticlesFluid::createInstance(numParticles, HOST);
	float* pxs = particles->m_pos->xs(HOST);
	float* pys = particles->m_pos->ys(HOST);
	float* pzs = particles->m_pos->zs(HOST);

	unsigned int id = 0;
	for (int i = -kernelRadiusPerEquilibriumDistance; i <= kernelRadiusPerEquilibriumDistance; ++i)
	{
		for (int j = -kernelRadiusPerEquilibriumDistance; j <= kernelRadiusPerEquilibriumDistance; ++j)
		{
			for (int k = -kernelRadiusPerEquilibriumDistance; k <= kernelRadiusPerEquilibriumDistance; ++k)
			{
				assert((id != centerPerticleId) || (id == centerPerticleId) && i == 0 && j== 0 && k == 0 );
				pxs[id] = (float)i * equilibriumDistance;
				pys[id] = (float)j * equilibriumDistance;
				pzs[id] = (float)k * equilibriumDistance;
				++id;
			}
		}
	}

	//Compute dotSumGradW and sumDotGradW.
	float sumDotGradW = 0.0f;
	for (unsigned int id = 0; id < numParticles; ++id)
	{
		float gradW[3];
		m_sphKernel.gradW(gradW, 0.0f, 0.0f, 0.0f, pxs[id], pys[id], pzs[id]);
		sumDotGradW += gradW[0] * gradW[0] + gradW[1] * gradW[1] + gradW[2] * gradW[2];
	}

	//Compute m_delta.
	float beta = m_deltaT * m_deltaT * m_particleMass * m_particleMass * 2.0f / (Constants::RO0 * Constants::RO0);
	float denominator = beta * sumDotGradW;
	m_delta = 1.0f / denominator;

	delete particles;
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void CalculatorPressurePciSph::calculation_host_(ParticlesFluid& particles, const CellCodeCalculator& ccc, const CompactHash& cHash)
{
	particles.sync(HOST);
	unsigned int size = particles.size();

	float* axs = particles.m_acceleration->xs(HOST);
	float* ays = particles.m_acceleration->ys(HOST);
	float* azs = particles.m_acceleration->zs(HOST);

	const float* vxs = particles.m_velocity->xs(HOST);
	const float* vys = particles.m_velocity->ys(HOST);
	const float* vzs = particles.m_velocity->zs(HOST);

	const float* pxs = particles.m_pos->xs(HOST);
	const float* pys = particles.m_pos->ys(HOST);
	const float* pzs = particles.m_pos->zs(HOST);

	//Temporary buffers for prediction-correction.
	PointSet predictedPosSet("predictedPosSet", size);
	predictedPosSet.allocate(HOST);
	float* ppxs = predictedPosSet.xs(HOST);
	float* ppys = predictedPosSet.ys(HOST);
	float* ppzs = predictedPosSet.zs(HOST);

	BufferFloat predictedPressureSet("predictedPressureSet", size);
	predictedPressureSet.allocate(HOST);
	predictedPressureSet.memset(0, HOST);
	float* pps = predictedPressureSet.get(HOST);

	BufferFloat predictedDensitySet("predictedDensitySet", size);
	predictedDensitySet.allocate(HOST);
	float* pds = predictedDensitySet.get(HOST);

	const unsigned int* sortedIdMaps = particles.m_sortedIdMap->get(HOST);

	float maxDensityError = FLT_MAX;
	for (unsigned int iter = 0; iter < m_numMaxIterations && maxDensityError / Constants::RO0 > m_maxRelativeDensityError; ++iter)
	{
		//Compute predicted fluid particle positions (time integral).
		for (unsigned int idP = 0; idP < size; ++idP)
		{
			//Implicit Euler.
			Point deltaXYZ;
			deltaXYZ[0] = (vxs[idP] + axs[idP] * m_deltaT) * m_deltaT;
			deltaXYZ[1] = (vys[idP] + ays[idP] * m_deltaT) * m_deltaT;
			deltaXYZ[2] = (vzs[idP] + azs[idP] * m_deltaT) * m_deltaT;

			ppxs[idP] = pxs[idP] + deltaXYZ[0] * m_deltaT;
			ppys[idP] = pys[idP] + deltaXYZ[1] * m_deltaT;
			ppzs[idP] = pzs[idP] + deltaXYZ[2] * m_deltaT;
		}

		//Compute density, density error -> update pressure to cancel the density error.
		maxDensityError = 0.0f;
		for (unsigned int idP = 0; idP < size; ++idP)
		{
			float sumW = 0.0f;

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
					const float distx = ppxs[idN] - ppxs[idP];
					const float disty = ppys[idN] - ppys[idP];
					const float distz = ppzs[idN] - ppzs[idP];
					const float dist2 = distx * distx + disty * disty + distz * distz;
					sumW += m_sphKernel.w(dist2);
				}
			}

			pds[idP] = m_particleMass * sumW;
			float densityError = pds[idP] - Constants::RO0;

			//The fluid gets glow up without this due to particle deficiency near the free surface.
			//Let's keep it non-negative until I implement e.g. ghost SPH.
			if (densityError < 0.0f)
			{
				densityError = 0.0f;
			}

			maxDensityError = std::max(maxDensityError, fabs(densityError));
			pps[idP] += m_delta * densityError; //Pressure to be canceled.
		}

		//Compute acceleration caused by the pressure.
		for (unsigned int idP = 0; idP < size; ++idP)
		{
			Point sumGradW(0.0f, 0.0f, 0.0f);
			Point sumCGradW(0.0f, 0.0f, 0.0f);

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
					float gradW[3];
					m_sphKernel.gradW(gradW, ppxs[idP], ppys[idP], ppzs[idP], ppxs[idN], ppys[idN], ppzs[idN]);
					sumGradW.x() += gradW[0];
					sumGradW.y() += gradW[1];
					sumGradW.z() += gradW[2];
					float c = pps[idN] / (pds[idN] * pds[idN]);
					sumCGradW.x() += c * gradW[0];
					sumCGradW.y() += c * gradW[1];
					sumCGradW.z() += c * gradW[2];
				}
			}

			const float c = pps[idP] / (pds[idP] * pds[idP]);


			Point currentCorrenctedAcceleration = - m_particleMass * (c * sumGradW + sumCGradW);

			//Acceleration limit. It prevents near free-surface partcles from moving too fast and corrupt the simulation.
			//Got the idea from OpenWorm (http://www.openworm.org/).
			//I wonder if there's a way to get rid of this.
			static const float MAX_DELTA_POS = 0.3f;
			const float maxAcceleration = MAX_DELTA_POS / (m_deltaT * m_deltaT); //Derived from PCISPH paper eq. (3).

			if (currentCorrenctedAcceleration.squaredNorm() > maxAcceleration * maxAcceleration)
			{
				currentCorrenctedAcceleration /= currentCorrenctedAcceleration.norm();
				currentCorrenctedAcceleration*= maxAcceleration;
			}

			axs[idP] += currentCorrenctedAcceleration.x();
			ays[idP] += currentCorrenctedAcceleration.y();
			azs[idP] += currentCorrenctedAcceleration.z();

		}
	}

	particles.setClean(HOST);
}
