#ifndef hohe_ParticlesFluid_H
#define hohe_ParticlesFluid_H

#include "ParticlesSph.h"

namespace hohehohe2
{

//-------------------------------------------------------------------
//-------------------------------------------------------------------
//Fluid Particles for SPH solver.
struct ParticlesFluid : public ParticlesSph
{

	///Particles density.
	BufferFloat* m_density;

	///Constructor.
	ParticlesFluid(const std::string& name="ParticlesFluid") : ParticlesSph(name), m_density(NULL){}

	///Destructor.
	virtual ~ParticlesFluid(){}

	virtual unsigned int size() const {return (m_pos)? m_pos->size() : 0;}

	///Create and setup for simulation. m_pos and m_velocity must be filled with initial values before the first step().
	static ParticlesFluid* createInstance(unsigned int size=0, MemoryType allocMemoryType=HOST)
	{
		ParticlesFluid* obj = new ParticlesFluid;
		fillMembers_(obj, size, allocMemoryType);
		return obj;
	};

protected:

	static void fillMembers_(ParticlesFluid* obj, unsigned int size, MemoryType allocMemoryType)
	{
		ParticlesSph::fillMembers_(obj, size, allocMemoryType);
		obj->addChild(obj->m_density = new BufferFloat("particle density"));
		obj->setSize(size);
		obj->allocate(allocMemoryType);
	}

};

}

#endif
