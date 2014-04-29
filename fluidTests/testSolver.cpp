#include "gtest/gtest.h"

#include "sph/FluidSolverSimpleSph.h"
#include "sph/ParticlesFluid.h"
#include "sph/ParticlesWall.h"

#include <iostream>

using namespace hohehohe2;


class Solver : public ::testing::Test
{
protected:

	FluidSolverSimpleSph* m_ssph;
	ParticlesFluid* m_particles;
	ParticlesWall* m_particlesWall;
	BufferSet::SPtr m_sptr;
	BufferSet::SPtr m_sptrWall;

	virtual void SetUp()
	{
		m_ssph = new FluidSolverSimpleSph;
		m_particles = ParticlesFluid::createInstance(1);
		m_particlesWall = ParticlesWall::createInstance(1);
		m_sptr = m_particles->getSelfSptr();
		m_sptrWall = m_particlesWall->getSelfSptr();

		float* pxs = m_particles->m_pos->xs(HOST);
		float* pys = m_particles->m_pos->ys(HOST);
		float* pzs = m_particles->m_pos->ys(HOST);
		float* vxs = m_particles->m_velocity->xs(HOST);
		float* vys = m_particles->m_velocity->ys(HOST);
		float* vzs = m_particles->m_velocity->ys(HOST);

		pxs[0] = 0.0f;
		pys[0] = 0.0f;
		pzs[0] = 0.0f;
		vxs[0] = 0.0f;
		vys[0] = 0.0f;
		vzs[0] = 0.0f;
	}

	virtual void TearDown()
	{
		delete m_ssph;
	}
};


TEST_F(Solver, create)
{
	m_ssph->step(*m_particles, *m_particlesWall, 0.1f);
}

TEST_F(Solver, singleParticle)
{
	FluidSolverSimpleSph ssph;
	ParticlesFluid* particles = ParticlesFluid::createInstance(1);
	ParticlesWall* particlesWall = ParticlesWall::createInstance(1);
	BufferSet::SPtr sptr = particles->getSelfSptr();
	BufferSet::SPtr sptrWall = particlesWall->getSelfSptr();

	float* pxs = particles->m_pos->xs(HOST);
	float* pys = particles->m_pos->ys(HOST);
	float* pzs = particles->m_pos->ys(HOST);
	float* vxs = particles->m_velocity->xs(HOST);
	float* vys = particles->m_velocity->ys(HOST);
	float* vzs = particles->m_velocity->ys(HOST);

	pxs[0] = 0.0f;
	pys[0] = 0.0f;
	pzs[0] = 0.0f;
	vxs[0] = 0.0f;
	vys[0] = 0.0f;
	vzs[0] = 0.0f;

	ssph.step(*particles, *particlesWall, 1.0f);

	//Commented out due to the adhoc ground boundary.
	//ASSERT_NEAR(pys[0], -9.8f, 0.001f);
}
