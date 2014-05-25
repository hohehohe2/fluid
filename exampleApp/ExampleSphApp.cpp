#include "ExampleSphApp.h"

#include <GL/glut.h>
//#include <hohe2Common/gl/Drawer.h>
#include "InitParticleDistributor.h"
#include <sph/ParticlesFluid.h>
#include <sph/ParticlesWall.h>

using namespace hohehohe2;


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void ExampleSphApp::reset(unsigned int id)
{
	PointSet* pos = new PointSet("particles pos");
	PointSet* velocity = new PointSet("particles velocity");
	PointSet* posWall = new PointSet("particlesWall pos");

	InitParticleDistributor::set(*pos, *velocity, *posWall, m_ssph.equilibriumDistance(), id);

	m_particles = ParticlesFluid::createInstance(true, pos->size());
	m_particles->setPos(pos);
	m_particles->setVelocity(velocity);
    m_sptr = m_particles->getSelfSptr();

	m_particlesWall = ParticlesWall::createInstance(true, posWall->size());
	m_particlesWall->setPos(posWall);
	m_sptrWall = m_particlesWall->getSelfSptr();
}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void ExampleSphApp::onKey(unsigned char key)
{
	switch (key)
	{
	case 27:
		exit(1);
	case 'r':
	case 'R':
		m_isSphereDraw = ! m_isSphereDraw;
		return;
	case ' ':
        m_ssph.step(*m_particles, *m_particlesWall, 0.01f);
		return;
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
	const float MAX_DENSITY = 1500.0f;
	const Point FLUID_RGB(1.0f, 0.0f, 0.0f);
	const Point WALL_RGB(0.5f, 0.5f, 0.5f);

	glPointSize(5.0f);

	float ed = m_ssph.equilibriumDistance();

	//----Fluid.
	m_particles->m_pos->sync(HOST);
	m_particles->m_density->sync(HOST);

	const float* pxs = m_particles->m_pos->xs(HOST);
	const float* pys = m_particles->m_pos->ys(HOST);
	const float* pzs = m_particles->m_pos->zs(HOST);
	const float* ds = m_particles->m_density->get(HOST);

	//glBegin(GL_POINTS);
	for (unsigned int i = 0; i < m_particles->m_pos->size(); ++i)
	{
		float scale = ds[i] / MAX_DENSITY;
		glColor3f(FLUID_RGB.x() * scale, FLUID_RGB.y() * scale, FLUID_RGB.z() * scale);
		if (i == 0)
		{
			glColor3f(0, 1, 1);
		}
		//glVertex3f(pxs[i], pys[i], pzs[i]);
		if (m_isSphereDraw)
		{
			glPushMatrix();
			glTranslatef(pxs[i], pys[i], pzs[i]);
			glutSolidSphere(ed / 2.0f, 16, 16);
			glPopMatrix();
		}
		else
		{
			glBegin(GL_POINTS);
			glVertex3f(pxs[i], pys[i], pzs[i]);
			glEnd();
		}
	}
	//glEnd();

	//----Wall.
	m_particlesWall->m_pos->sync(HOST);

	const float* wpxs = m_particlesWall->m_pos->xs(HOST);
	const float* wpys = m_particlesWall->m_pos->ys(HOST);
	const float* wpzs = m_particlesWall->m_pos->zs(HOST);

	glBegin(GL_POINTS);
	for (unsigned int i = 0; i < m_particlesWall->m_pos->size(); ++i)
	{
		float scale = 1.0f;
		glColor3f(WALL_RGB.x() * scale, WALL_RGB.y() * scale, WALL_RGB.z() * scale);
		glVertex3f(wpxs[i], wpys[i], wpzs[i]);
	}
	glEnd();
}
