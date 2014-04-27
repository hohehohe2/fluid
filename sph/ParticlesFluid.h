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

	///Particles velocity.
	PointSet* m_velocity;

	///Particles density.
	BufferFloat* m_density;

	///Constructor.
	ParticlesFluid(const std::string& name="ParticlesFluid") : ParticlesSph(name), m_velocity(NULL), m_density(NULL){}

	///Destructor.
	virtual ~ParticlesFluid(){}

	virtual unsigned int size() const {return (m_pos)? m_pos->size() : 0;}

	//Make sure every size is the same.
	void setVelocity(PointSet* velocity){removeChild(m_velocity); addChild(m_velocity = velocity);}
	void setAcceleration(PointSet* acceleration){removeChild(m_acceleration); addChild(m_acceleration = acceleration);}

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
		obj->addChild(obj->m_velocity = new PointSet("particle velocity"));
		obj->addChild(obj->m_density = new BufferFloat("particle density"));
		obj->setSize(size);
		obj->allocate(allocMemoryType);
		obj->m_velocity->memset(0, HOST);
	}

};

}

#endif
