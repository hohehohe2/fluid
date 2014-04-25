#include "ExampleSphApp.h"
#include <hohe2Common/gl/Drawer.h>

using namespace hohehohe2;


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void ExampleSphApp::reset(unsigned int type)
{
    if (type == 0)
    {
        m_particles = FluidSolverSimpleSph::Particles::createInstance(1);
        m_sptr = m_particles->getSelfSptr();

		float* pxs = m_particles->m_pos->xs(true);
		float* pys = m_particles->m_pos->ys(true);
		float* vxs = m_particles->m_velocity->xs(true);
		float* vys = m_particles->m_velocity->ys(true);

		pxs[0] = 0.0f;
		pys[0] = 0.0f;
		vxs[0] = 0.0f;
		vys[0] = 0.0f;
	}
}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void ExampleSphApp::onKey(unsigned char key)
{
    if (key == ' ')
    {
        m_ssph.step(*m_particles, 0.1f);
    }
    else if (key == '0')
    {
        reset(0);
    }
}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void ExampleSphApp::draw()
{
	Drawer::draw(*m_particles->m_pos);
}
