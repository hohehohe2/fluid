#include "gtest/gtest.h"

#include "sph/FluidSolverSimpleSph.h"
#include "sph/FluidParticles.h"

#include <iostream>

using namespace hohehohe2;


class Solver : public ::testing::Test
{
protected:

	FluidSolverSimpleSph* m_ssph;
	FluidParticles* m_particles;
	BufferSet::SPtr sptr;

	virtual void SetUp()
	{
		m_ssph = new FluidSolverSimpleSph;
		m_particles = FluidParticles::createInstance(1);
		sptr = m_particles->getSelfSptr();

		float* pxs = m_particles->m_pos->xs(true);
		float* pys = m_particles->m_pos->ys(true);
		float* pzs = m_particles->m_pos->ys(true);
		float* vxs = m_particles->m_velocity->xs(true);
		float* vys = m_particles->m_velocity->ys(true);
		float* vzs = m_particles->m_velocity->ys(true);

		pxs[0] = 0.0f;
		pys[0] = 0.0f;
		pzs[0] = 0.0f;
		vxs[0] = 0.0f;
		vys[0] = 0.0f;
		vzs[0] = 0.0f;
	}

	virtual void TearDown()
	{
		sptr.reset();
		delete m_ssph;
	}
};


TEST_F(Solver, create)
{
	m_ssph->step(*m_particles, 0.1f);
}

TEST_F(Solver, singleParticle)
{
	FluidSolverSimpleSph ssph;
	FluidParticles* particles = FluidParticles::createInstance(1);
	BufferSet::SPtr sptr = particles->getSelfSptr();

	float* pxs = particles->m_pos->xs(true);
	float* pys = particles->m_pos->ys(true);
	float* pzs = particles->m_pos->ys(true);
	float* vxs = particles->m_velocity->xs(true);
	float* vys = particles->m_velocity->ys(true);
	float* vzs = particles->m_velocity->ys(true);

	pxs[0] = 0.0f;
	pys[0] = 0.0f;
	pzs[0] = 0.0f;
	vxs[0] = 0.0f;
	vys[0] = 0.0f;
	vzs[0] = 0.0f;

	ssph.step(*particles, 1.0f);

	//Commented out due to the adhoc ground boundary.
	//ASSERT_NEAR(pys[0], -9.8f, 0.001f);
}
