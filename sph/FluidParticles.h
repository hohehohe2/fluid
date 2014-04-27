#ifndef hohe_FluidParticles_H
#define hohe_FluidParticles_H

#include <hohe2Common/geo/basicGeos.h>

namespace hohehohe2
{

//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Particles for SPH solver.
struct FluidParticles : public BufferSetSized
{

	///Particles position. Always exists.
	PointSet* m_pos;

	///Particles velocity.
	PointSet* m_velocity;

	///Particles acceleration.
	PointSet* m_acceleration;

	///Particles density.
	BufferFloat* m_density;

	///m_pos etc. index sorted by its code.
	BufferUInt* m_sortedIdMap;

	///Constructor.
	FluidParticles(const std::string& name="FluidParticles") : BufferSetSized(name), m_pos(NULL), m_velocity(NULL), m_acceleration(NULL), m_density(NULL), m_sortedIdMap(NULL){}

	///Destructor.
	virtual ~FluidParticles(){}

	virtual unsigned int size() const {return (m_pos)? m_pos->size() : 0;}

	//Make sure every size is the same.
	void setPos(PointSet* pos){removeChild(m_pos); addChild(m_pos = pos);}
	void setVelocity(PointSet* velocity){removeChild(m_velocity); addChild(m_velocity = velocity);}
	void setAcceleration(PointSet* acceleration){removeChild(m_acceleration); addChild(m_acceleration = acceleration);}
	void setDensity(BufferFloat* density){removeChild(m_density); addChild(m_density = density);}
	void setSortedIdMap(BufferUInt* sortedIdMap){removeChild(m_sortedIdMap); addChild(m_sortedIdMap = sortedIdMap);}

	///Create and setup for simulation. m_pos and m_velocity must be filled with initial values before the first step().
	static FluidParticles* createInstance(unsigned int size=0, MemoryType allocMemoryType=HOST)
	{
		//Separated object creation and filling member variables so that we can reuse fillMembers_() in a subclass.
		FluidParticles* obj = new FluidParticles;
		fillMembers_(obj, size, allocMemoryType);
		return obj;
	};

protected:

	static void fillMembers_(FluidParticles* obj, unsigned int size, MemoryType allocMemoryType)
	{
		obj->addChild(obj->m_pos = new PointSet("particle pos"));
		obj->addChild(obj->m_velocity = new PointSet("particle velocity"));
		obj->addChild(obj->m_acceleration = new PointSet("particle acceleration"));
		obj->addChild(obj->m_density = new BufferFloat("particle density"));
		obj->addChild(obj->m_sortedIdMap = new BufferUInt("particle sortedIdMap"));
		obj->setSize(size);
		obj->allocate(allocMemoryType);
		obj->m_pos->memset(0, HOST);
		obj->m_velocity->memset(0, HOST);
	}

};

}

#endif
