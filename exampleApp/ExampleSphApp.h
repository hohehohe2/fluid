#ifndef hohe_ExampleSphApp_H
#define hohe_ExampleSphApp_H

#include <sph/FluidSolverSimpleSph.h>

namespace hohehohe2
{

class ExampleSphApp
{
public:

	ExampleSphApp() : m_isSphereDraw(true){}

	void reset(unsigned int id=0);
	void draw();
	void onKey(unsigned char key);

private:

	FluidSolverSimpleSph m_ssph;
	ParticlesFluid* m_particles;
	ParticlesWall* m_particlesWall;
	BufferSet::SPtr m_sptr;
	BufferSet::SPtr m_sptrWall;

	bool m_isSphereDraw;
};

}

#endif
