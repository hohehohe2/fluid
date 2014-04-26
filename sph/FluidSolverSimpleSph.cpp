#include "FluidSolverSimpleSph.h"
#include "Constants.h"

using namespace hohehohe2;


//-------------------------------------------------------------------
//-------------------------------------------------------------------
const float FluidSolverSimpleSph::K = 20000.0f;
const float FluidSolverSimpleSph::PET_PEEVE_COURANT_NUMBER = 0.1f;


//-------------------------------------------------------------------
//-------------------------------------------------------------------
FluidSolverSimpleSph::FluidSolverSimpleSph(float particleMass)
	:
	m_particleMass(particleMass),
	m_cHash(COMPACT_HASH_NUM_HASH_ENTRIES, COMPACT_HASH_NUM_ELEMENTS_IN_A_LIST, COMPACT_HASH_NUM_LISTS, HOST)
{
	//Adjust the kernel radius so that several dozens of neighbor particles are in the radius at rest density.
	float particleVolume = particleMass / Constants::RO0;
	m_restLength = (float)pow(particleVolume, 1.0 / 3.0);
	double kernelRadius = m_restLength * 4.0;
	m_sphKernel.setR((float)kernelRadius);
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void FluidSolverSimpleSph::step(Particles& particles, float deltaT)
{
	float remaining = deltaT;
	bool loop = true;
	do
	{
		float maxVelocity = calcMaxVelocity_(particles);
		float dt = PET_PEEVE_COURANT_NUMBER / maxVelocity;
		if (dt > remaining)
		{
			dt = remaining;
			loop = false;
		}
		remaining -= dt;

		updateNeighbors_(particles);
		calcDensity_(particles);
		calcAcceleration_(particles);
		integrate_(particles, dt);

	} while(loop);
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
float FluidSolverSimpleSph::calcMaxVelocity_(Particles& particles)
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
void FluidSolverSimpleSph::updateNeighbors_(Particles& particles)
{
	particles.sync(HOST);
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void FluidSolverSimpleSph::calcDensity_(Particles& particles)
{
	particles.sync(HOST);

	const float r2 = m_sphKernel.r() * m_sphKernel.r();

	particles.m_pos->sync(HOST);
	const float* pxs = particles.m_pos->xs(HOST);
	const float* pys = particles.m_pos->ys(HOST);
	const float* pzs = particles.m_pos->zs(HOST);
	float* ds = particles.m_density->get(HOST);

	unsigned int size = particles.size();

	for (unsigned int idP = 0; idP < size; ++idP)
	{
		ds[idP] = 0.0f;
		for (unsigned int idN = 0; idN < size; ++idN)
		{
			float distx = pxs[idN] - pxs[idP];
			float disty = pys[idN] - pys[idP];
			float distz = pzs[idN] - pzs[idP];
			float dist2 = distx * distx + disty * disty + distz * distz;
			ds[idP] += m_sphKernel.w(dist2);
		}
		ds[idP] *= m_particleMass;

		//Set rest density as minimum to fake air pressure.
		if (ds[idP] < Constants::RO0)
		{
			ds[idP] = Constants::RO0;
		}
	}
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void FluidSolverSimpleSph::calcAcceleration_(Particles& particles)
{
	particles.sync(HOST);

	particles.m_pos->sync(HOST);
	particles.m_density->sync(HOST);
	const float* pxs = particles.m_pos->xs(HOST);
	const float* pys = particles.m_pos->ys(HOST);
	const float* pzs = particles.m_pos->zs(HOST);
	const float* ds = particles.m_density->get(HOST);
	float* axs = particles.m_acceleration->xs(HOST);
	float* ays = particles.m_acceleration->ys(HOST);
	float* azs = particles.m_acceleration->zs(HOST);

	unsigned int size = particles.size();

	for (unsigned int idP = 0; idP < size; ++idP)
	{
		float densityP = ds[idP];

		///Gravity.
		axs[idP] = 0.0f;
		ays[idP] = Constants::G;
		azs[idP] = 0.0f;

		///Pressure.
		for (unsigned int idN = 0; idN < size; ++idN)
		{
			float gradW[3];
			m_sphKernel.gradW(gradW, pxs[idP], pys[idP], pzs[idP], pxs[idN], pys[idN], pzs[idN]);

			float densityN = ds[idN];
			float pressureP = densityToPressure_(densityP);
			float pressureN = densityToPressure_(densityN);
			float c = m_particleMass * (pressureP + pressureN) / (2.0f * densityN);
			float gradPx = c * gradW[0];
			float gradPy = c * gradW[1];
			float gradPz = c * gradW[2];
			axs[idP] += - gradPx / densityP;
			ays[idP] += - gradPy / densityP;
			azs[idP] += - gradPz / densityP;
		}

		///Viscosity omitted.
	}
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void FluidSolverSimpleSph::integrate_(Particles& particles, float deltaT)
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
	}
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
float FluidSolverSimpleSph::densityToPressure_(float density)
{
	return K * (density * density / (Constants::RO0 * Constants::RO0) - 1.0f);

}
