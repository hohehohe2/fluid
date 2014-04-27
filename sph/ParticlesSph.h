#ifndef hohe_ParticlesSph_H
#define hohe_ParticlesSph_H

#include <hohe2Common/geo/basicGeos.h>

namespace hohehohe2
{

//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Particles for SPH solver.
struct ParticlesSph : public BufferSetSized
{

	///Particles position. Always exists.
	PointSet* m_pos;

	///Particles acceleration.
	PointSet* m_acceleration;

	///m_pos etc. index sorted by its code.
	BufferUInt* m_sortedIdMap;

	///Constructor.
	ParticlesSph(const std::string& name="ParticlesSph") : BufferSetSized(name), m_pos(NULL), m_acceleration(NULL), m_sortedIdMap(NULL){}

	///Destructor.
	virtual ~ParticlesSph(){}

	virtual unsigned int size() const {return (m_pos)? m_pos->size() : 0;}

	//Make sure every size is the same.
	void setPos(PointSet* pos){removeChild(m_pos); addChild(m_pos = pos);}
	void setAcceleration(PointSet* acceleration){removeChild(m_acceleration); addChild(m_acceleration = acceleration);}
	void setSortedIdMap(BufferUInt* sortedIdMap){removeChild(m_sortedIdMap); addChild(m_sortedIdMap = sortedIdMap);}

	///Create and setup for simulation. m_pos and m_velocity must be filled with initial values before the first step().
	static ParticlesSph* createInstance(unsigned int size=0, MemoryType allocMemoryType=HOST)
	{
		ParticlesSph* obj = new ParticlesSph;
		fillMembers_(obj, size, allocMemoryType);
		return obj;
	};

protected:

	static void fillMembers_(ParticlesSph* obj, unsigned int size, MemoryType allocMemoryType)
	{
		obj->addChild(obj->m_pos = new PointSet("particle pos"));
		obj->addChild(obj->m_acceleration = new PointSet("particle acceleration"));
		obj->addChild(obj->m_sortedIdMap = new BufferUInt("particle sortedIdMap"));
		obj->setSize(size);
		obj->allocate(allocMemoryType);
		obj->m_pos->memset(0, HOST);
	}

};

}

#endif
