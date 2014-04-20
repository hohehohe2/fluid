#ifndef hohe_Particles_H
#define hohe_Particles_H

#include <cudaCommon/Buffer.h>
#include "basicGeos.h"

namespace hohehohe2
{

//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Particle cloud.
class Particles : public BufferSetSameSizedChildren
{

public:

	typedef std::shared_ptr < Particles > SPtr;

	///Constructor.
	Particles(const std::string& name="", unsigned int size=0, MemoryType allocMemoryType=BufferSet::HOST) : BufferSetSameSizedChildren(name)
	{
		addChild(m_ppPos = new Points("particle ppPos"));
		addChild(m_pPos = new Points("particle pPos"));
		addChild(m_pos = new Points("particle pos"));
		addChild(m_velocity = new BufferFloat("particle velocity"));
		addChild(m_mass = new BufferFloat("particle mass"));
		addChild(m_force = new Points("particle force"));
		addChild(m_pressure = new BufferFloat("particle pressure"));
		addChild(m_density = new BufferFloat("particle density"));
		setSize(size);
		allocate(allocMemoryType);
	}

	///Destructor.
	virtual ~Particles(){}

	virtual unsigned int size() const {return m_pos->size();}

public:

	///Particle prev prev position.
	Points* m_ppPos;

	///Particle prev position.
	Points* m_pPos;

	///Particle position. Always exists.
	Points* m_pos;

	///Particle velocity.
	BufferFloat* m_velocity;

	///Particle mass.
	BufferFloat* m_mass;

	///Particle force.
	Points* m_force;

	///Particle pressure.
	BufferFloat* m_pressure;

	///Particle density.
	BufferFloat* m_density;

};

}

#endif
