#include "FluidSolverSimpleSph.h"
#include <hohe2Common/util/BufferUtil.h>
#include <hohe2Common/container/CellCodeCalculator.h>
#include "Constants.h"
#include "FluidParticles.h"
#include "PressureCalculator.h"


using namespace hohehohe2;


//-------------------------------------------------------------------
//-------------------------------------------------------------------
const float FluidSolverSimpleSph::MU = 0.1f;
const float FluidSolverSimpleSph::PET_PEEVE_COURANT_NUMBER = 0.5f;

//-------------------------------------------------------------------
//-------------------------------------------------------------------
FluidSolverSimpleSph::FluidSolverSimpleSph(float particleMass)
	:
	m_particleMass(particleMass),
	m_cHash(COMPACT_HASH_NUM_HASH_ENTRIES, COMPACT_HASH_NUM_ELEMENTS_IN_A_LIST, COMPACT_HASH_NUM_LISTS, HOST),
	m_pressureCalculator(particleMass)
{
	//Adjust the kernel radius so that several dozens of neighbor particles are in the radius at rest density.
	float particleVolume = particleMass / Constants::RO0;
	m_restLength = (float)pow(particleVolume, 1.0 / 3.0);
	double kernelRadius = m_restLength * 4.0;
	m_sphKernel.setR((float)kernelRadius);
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void FluidSolverSimpleSph::step(FluidParticles& particles, float deltaT)
{
	float remaining = deltaT;
	bool loop = true;
	do
	{
		float maxVelocity = calcMaxVelocity_(particles);
		float dt = PET_PEEVE_COURANT_NUMBER * m_sphKernel.r() / maxVelocity;
		if (dt > remaining)
		{
			dt = remaining;
			loop = false;
		}
		remaining -= dt;

		std::cout << remaining << ": updateNeighbors - ";
		updateNeighbors_(particles);
		std::cout << "calcDensity - ";
		calcDensity_host_(particles);
		std::cout << "calcAcceleration - ";
		calcAcceleration_host_(particles);
		m_pressureCalculator.calcAcceleration(particles, m_sphKernel, HOST, m_cHash);
		std::cout << "integrate\n";
		integrate_(particles, dt);

	} while(loop);
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
float FluidSolverSimpleSph::calcMaxVelocity_(const FluidParticles& particles)
{
	particles.sync(HOST);

	float maxVelocity2 = 0;

	particles.m_velocity->sync(HOST);
	const float* vxs = particles.m_velocity->xs(HOST);
	const float* vys = particles.m_velocity->ys(HOST);
	const float* vzs = particles.m_velocity->zs(HOST);

	unsigned int size = particles.size();
	for (unsigned int idP = 0; idP < size; ++idP)
	{
		float mv2 = vxs[idP] * vxs[idP] + vys[idP] * vys[idP] + vzs[idP] * vzs[idP];
		if (maxVelocity2 < mv2)
		{
			maxVelocity2 = mv2;
		}
	}
	return sqrt(maxVelocity2);
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void FluidSolverSimpleSph::updateNeighbors_(FluidParticles& particles)
{
	const float kernelRaidus = m_sphKernel.r();

	particles.sync(HOST);

	BoundingBox bbox;
	particles.m_pos->calcBoundingBox(bbox, HOST);

	//Shift min values a little bit to make sure the particles which have min values are inside the cell.
	bbox.m_min -= Point(kernelRaidus, kernelRaidus, kernelRaidus) / 100.0f;
	particles.m_pos->m_lastCalculatedBbox.m_min = bbox.m_min;

	CellCodeCalculator ccc(bbox.m_min, kernelRaidus);

	//Get the unsorted code list.
	BufferUInt codeSet;
	ccc.getCode32(codeSet, *particles.m_pos, HOST);

	//Sort it. Once we create the hash, we don't have to keep using it so it's not a particles member variable.
	BufferUInt sortedCodeSet;
	BufferUtil::sortByKey(*particles.m_sortedIdMap, sortedCodeSet, codeSet, DEVICE);

	//Create compact hash.
	if ( ! m_cHash.build(sortedCodeSet, HOST))
	{
		std::cerr << "Hash full\n";
	}
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void FluidSolverSimpleSph::calcDensity_host_(FluidParticles& particles)
{
	particles.sync(HOST);

	particles.m_pos->sync(HOST);
	particles.m_sortedIdMap->sync(HOST);
	const float* pxs = particles.m_pos->xs(HOST);
	const float* pys = particles.m_pos->ys(HOST);
	const float* pzs = particles.m_pos->zs(HOST);
	const unsigned int* sortedIdMaps = particles.m_sortedIdMap->get(HOST);

	float* ds = particles.m_density->get(HOST);

	unsigned int size = particles.size();

	CellCodeCalculator ccc(particles.m_pos->m_lastCalculatedBbox.m_min, m_sphKernel.r());

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
			unsigned int numObjects = m_cHash.lookup(index, code);
			for (unsigned int j = 0; j < numObjects; ++j)
			{
				unsigned int idN = sortedIdMaps[index + j];
				float distx = pxs[idN] - pxs[idP];
				float disty = pys[idN] - pys[idP];
				float distz = pzs[idN] - pzs[idP];
				float dist2 = distx * distx + disty * disty + distz * distz;
				ds[idP] += m_sphKernel.w(dist2);
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


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void FluidSolverSimpleSph::calcAcceleration_host_(FluidParticles& particles)
{
	particles.sync(HOST);

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

	unsigned int size = particles.size();

	CellCodeCalculator ccc(particles.m_pos->m_lastCalculatedBbox.m_min, m_sphKernel.r());

	#pragma omp parallel for
	for (int idP = 0; idP < (int)size; ++idP)
	{

		//----Gravity.
		axs[idP] = 0.0f;
		ays[idP] = Constants::G;
		azs[idP] = 0.0f;

		//----Viscosity.
		for (unsigned int i= 0; i < 27; ++i)
		{
			bool isValid;
			unsigned int code = ccc.getNeighborCode32(isValid, pxs[idP], pys[idP], pzs[idP], i);
			if ( ! isValid)
			{
				continue;
			}

			unsigned int index;
			unsigned int numObjects = m_cHash.lookup(index, code);
			for (unsigned int j = 0; j < numObjects; ++j)
			{
				unsigned int idN = sortedIdMaps[index + j];

				float distx = pxs[idN] - pxs[idP];
				float disty = pys[idN] - pys[idP];
				float distz = pzs[idN] - pzs[idP];
				float dist2 = distx * distx + disty * disty + distz * distz;
				float laplace = m_sphKernel.laplaceW(dist2);
				float viscosCoef = MU / ds[idN] * laplace;

                axs[idP] += viscosCoef * (vxs[idN] - vxs[idP]) / m_particleMass;
                ays[idP] += viscosCoef * (vys[idN] - vys[idP]) / m_particleMass;
                azs[idP] += viscosCoef * (vzs[idN] - vzs[idP]) / m_particleMass;
			}
		}

	}

	particles.setClean(HOST);
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void FluidSolverSimpleSph::integrate_(FluidParticles& particles, float deltaT)
{
	particles.sync(HOST);

	particles.m_acceleration->sync(HOST);
	const float* axs = particles.m_acceleration->xs(HOST);
	const float* ays = particles.m_acceleration->ys(HOST);
	const float* azs = particles.m_acceleration->zs(HOST);
	float* pxs = particles.m_pos->xs(HOST);
	float* pys = particles.m_pos->ys(HOST);
	float* pzs = particles.m_pos->zs(HOST);
	float* vxs = particles.m_velocity->xs(HOST);
	float* vys = particles.m_velocity->ys(HOST);
	float* vzs = particles.m_velocity->zs(HOST);

	unsigned int size = particles.size();

	for (unsigned int idP = 0; idP < size; ++idP)
	{
		vxs[idP] += axs[idP] * deltaT;
		vys[idP] += ays[idP] * deltaT;
		vzs[idP] += azs[idP] * deltaT;
		pxs[idP] += vxs[idP] * deltaT;
		pys[idP] += vys[idP] * deltaT;
		pzs[idP] += vzs[idP] * deltaT;

		//Adhoc boundary.
		if (pxs[idP] < 0.0f)
		{
			pxs[idP] = 0.0f;
			vxs[idP] = 0;
		}

		if (pzs[idP] < 0.0f)
		{
			pzs[idP] = 0.0f;
			vzs[idP] = 0;
		}

		if (pys[idP] < 0.0f)
		{
			pys[idP] = 0.0f;
			vys[idP] = 0;
		}
	}

	particles.setClean(HOST);
}

