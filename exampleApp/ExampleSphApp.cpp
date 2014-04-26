#include "ExampleSphApp.h"

#include <hohe2Common/gl/Drawer.h>
#include "InitParticleDistributor.h"

using namespace hohehohe2;


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void ExampleSphApp::reset(unsigned int id)
{
	PointSet* pos = new PointSet("particle pos");
	PointSet* velocity = new PointSet("particle velocity");

	InitParticleDistributor::set(*pos, *velocity, m_ssph.restLength(), id);

	m_particles = FluidSolverSimpleSph::Particles::createInstance(pos->size());
	m_particles->setPos(pos);
	m_particles->setVelocity(velocity);
    m_sptr = m_particles->getSelfSptr();

}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void ExampleSphApp::onKey(unsigned char key)
{
    if (key == 27)
    {
		exit(1);
	}

	if (key == ' ')
    {
        m_ssph.step(*m_particles, 0.1f);
		return;
    }

	switch (key)
	{
	case '0':
	case '1':
	case '2':
	case '3':
	case '4':
	case '5':
	case '6':
	case '7':
	case '8':
	case '9':
        reset(key - '0');
		break;
	default:
		break;
	}
}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void ExampleSphApp::draw()
{
	Drawer::draw(*m_particles->m_pos);
}
