#ifndef hohe_ParticlesWall_H
#define hohe_ParticlesWall_H

#include "ParticlesSph.h"

namespace hohehohe2
{

//-------------------------------------------------------------------
//-------------------------------------------------------------------
//Wall Particles for SPH solver.
struct ParticlesWall : public ParticlesSph
{

	///Particle volumes.
	BufferFloat* m_volume;

	///Constructor.
	ParticlesWall(const std::string& name="ParticlesWall") : ParticlesSph(name), m_volume(NULL){}

	///Destructor.
	virtual ~ParticlesWall(){}

	virtual unsigned int size() const {return (m_pos)? m_pos->size() : 0;}

	//Make sure every size is the same.
	void setAcceleration(PointSet* acceleration){removeChild(m_acceleration); addChild(m_acceleration = acceleration);}

	///Create and setup for simulation. m_pos and m_velocity must be filled with initial values before the first step().
	static ParticlesWall* createInstance(unsigned int size=0, MemoryType allocMemoryType=HOST)
	{
		ParticlesWall* obj = new ParticlesWall;
		fillMembers_(obj, size, allocMemoryType);
		return obj;
	};

protected:

	static void fillMembers_(ParticlesWall* obj, unsigned int size, MemoryType allocMemoryType)
	{
		ParticlesSph::fillMembers_(obj, size, allocMemoryType);
		obj->addChild(obj->m_volume = new BufferFloat("particle volume"));
		obj->setSize(size);
		obj->allocate(allocMemoryType);
	}

};

}

#endif
