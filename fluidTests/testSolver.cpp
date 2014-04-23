#include "gtest/gtest.h"

#include "sph/FluidSolverSimpleSph.h"

#include <iostream>

using namespace hohehohe2;


class Solver : public ::testing::Test
{
protected:

	FluidSolverSimpleSph* m_ssph;
	FluidSolverSimpleSph::Particles* particles;
	BufferSet::SPtr sptr;
 
    virtual void SetUp()
    {
		m_ssph = new FluidSolverSimpleSph;
		particles = FluidSolverSimpleSph::Particles::createInstance(1);
		sptr = particles->getSelfSptr();

		float* pxs = particles->m_pos->xs(true);
		float* pys = particles->m_pos->ys(true);
		float* vxs = particles->m_velocity->xs(true);
		float* vys = particles->m_velocity->ys(true);
    }

	virtual void TearDown()
    {
		sptr.reset();
        delete m_ssph;
    }
};


TEST_F(Solver, create)
{
	m_ssph->step(*particles, 0.1f);
}

TEST_F(Solver, singleParticle)
{
	FluidSolverSimpleSph ssph;
	FluidSolverSimpleSph::Particles* particles = FluidSolverSimpleSph::Particles::createInstance(1);
	BufferSet::SPtr sptr = particles->getSelfSptr();

	float* pxs = particles->m_pos->xs(true);
	float* pys = particles->m_pos->ys(true);
	float* vxs = particles->m_velocity->xs(true);
	float* vys = particles->m_velocity->ys(true);

	pxs[0] = 0.0f;
	pys[0] = 0.0f;
	vxs[0] = 0.0f;
	vys[0] = 0.0f;

	ssph.step(*particles, 1.0f);

	ASSERT_NEAR(pys[0], -9.8f, 0.001f);
}
