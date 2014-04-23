#include "FluidSolverSimpleSph.h"
#include "Constants.h"

using namespace hohehohe2;


//-------------------------------------------------------------------
//-------------------------------------------------------------------
const float FluidSolverSimpleSph::K = 20000.0f;
const float FluidSolverSimpleSph::PET_PEEVE_COURANT_NUMBER = 0.1f;


//-------------------------------------------------------------------
//-------------------------------------------------------------------
FluidSolverSimpleSph::FluidSolverSimpleSph(float particleMass) : m_particleMass(particleMass)
{
	//Adjust the kernel radius so that several dozens of neighbor particles are in the radius at rest density.
	float particleVolume = particleMass / Constants::RO0;
	double kernelRadius = pow(particleVolume, 1.0 / 3.0) * 4.0;
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
    particles.sync();

	float maxVelocity2 = 0;

	float* vxs = particles.m_velocity->xs(true);
	float* vys = particles.m_velocity->ys(true);
	float* vzs = particles.m_velocity->zs(true);

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
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void FluidSolverSimpleSph::calcDensity_(Particles& particles)
{
    particles.sync();

	const float r2 = m_sphKernel.r() * m_sphKernel.r();

	float* pxs = particles.m_pos->xs(true);
	float* pys = particles.m_pos->ys(true);
	float* pzs = particles.m_pos->zs(true);
	float* ds = particles.m_density->get(true);

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
    particles.sync();

	float* pxs = particles.m_pos->xs(true);
	float* pys = particles.m_pos->ys(true);
	float* pzs = particles.m_pos->zs(true);
	float* axs = particles.m_acceleration->xs(true);
	float* ays = particles.m_acceleration->ys(true);
	float* azs = particles.m_acceleration->zs(true);
	float* ds = particles.m_density->get(true);

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
    particles.sync();

	float* pxs = particles.m_pos->xs(true);
	float* pys = particles.m_pos->ys(true);
	float* pzs = particles.m_pos->zs(true);
	float* vxs = particles.m_velocity->xs(true);
	float* vys = particles.m_velocity->ys(true);
	float* vzs = particles.m_velocity->zs(true);
	float* axs = particles.m_acceleration->xs(true);
	float* ays = particles.m_acceleration->ys(true);
	float* azs = particles.m_acceleration->zs(true);

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
