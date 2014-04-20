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
	/**
	Only position is added to a child by default. Other children are added if specified as the parameters.
	**/
	Particles(const std::string& name="",
		bool ppPos=false, bool pPos=false, bool velocity=false, bool mass=false, bool force=false, bool pressure=false, bool density=false,
		unsigned int size=0, MemoryType allocMemoryType=BufferSet::HOST)
		: BufferSetSameSizedChildren(name)
	{
		createChildren(ppPos, pPos, velocity, mass, force, pressure, density, allocMemoryType);
	}

	///Destructor.
	virtual ~Particles(){}

	virtual unsigned int size() const {return m_pos->size();}

	///Create children.
	void createChildren(bool ppPos, bool pPos, bool velocity, bool mass, bool force, bool pressure, bool density,
		unsigned int size=0, MemoryType allocMemoryType=BufferSet::HOST)
	{
		if (ppPos)
		{
			removeChild(m_ppPos);
			addChild(m_ppPos = new Points("particle ppPos"));
		}
		if (pPos)
		{
			removeChild(m_pPos);
			addChild(m_pPos = new Points("particle pPos"));
		}
		if (velocity)
		{
			removeChild(m_velocity);
			addChild(m_velocity = new BufferFloat("particle velocity"));
		}
		if (mass)
		{
			removeChild(m_mass);
			addChild(m_mass = new BufferFloat("particle mass"));
		}
		if (force)
		{
			removeChild(m_force);
			addChild(m_force = new Points("particle force"));
		}
		if (pressure)
		{
			removeChild(m_pressure);
			addChild(m_pressure = new BufferFloat("particle pressure"));
		}
		if (density)
		{
			removeChild(m_density);
			addChild(m_density = new BufferFloat("particle density"));
		}
		setSize(size);
		allocate(allocMemoryType);
	}

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
