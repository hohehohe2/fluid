#include "CalculatorSemiImplicitEulerIntegrate.h"

#include "ParticlesFluid.h"

using namespace hohehohe2;


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void CalculatorSemiImplicitEulerIntegrate::integrate_host_(ParticlesFluid& particles, float deltaT)
{
	particles.m_acceleration->sync(HOST);
	particles.m_velocity->sync(HOST);
	particles.m_pos->sync(HOST);
	const float* axs = particles.m_acceleration->xs(HOST);
	const float* ays = particles.m_acceleration->ys(HOST);
	const float* azs = particles.m_acceleration->zs(HOST);
	float* vxs = particles.m_velocity->xs(HOST);
	float* vys = particles.m_velocity->ys(HOST);
	float* vzs = particles.m_velocity->zs(HOST);
	float* pxs = particles.m_pos->xs(HOST);
	float* pys = particles.m_pos->ys(HOST);
	float* pzs = particles.m_pos->zs(HOST);

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

	particles.setClean(HOST);
}
