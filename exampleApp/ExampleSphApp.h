#ifndef hohe_ExampleSphApp_H
#define hohe_ExampleSphApp_H

#include "sph/FluidSolverSimpleSph.h"

namespace hohehohe2
{

class ExampleSphApp
{
public:
	void reset(unsigned int type=0);
	void draw();
	void onKey(unsigned char key);

private:

	FluidSolverSimpleSph m_ssph;
	FluidSolverSimpleSph::Particles* m_particles;
	BufferSet::SPtr m_sptr;

};

}

#endif
