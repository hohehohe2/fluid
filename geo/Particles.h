#ifndef hohe_Particles_H
#define hohe_Particles_H

#include <cudaCommon/Buffer.h>
#include "basicGeos.h"

namespace hohehohe2
{

//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Particle cloud.
class Particles : public BufferSetSized
{

public:

	typedef std::shared_ptr < Particles > SPtr;

	///Constructor.
	/**
	Only position is added to a child by default. Other children are added if specified as the parameters.
	**/
	Particles(const std::string& name="",
		bool ppPos=false, bool pPos=false, bool velocity=false, bool acceleration=false, bool mass=false, bool force=false, bool pressure=false, bool density=false,
		unsigned int size=0, MemoryType allocMemoryType=BufferSet::HOST)
		: BufferSetSized(name)
	{
		createChildren(ppPos, pPos, velocity, acceleration, mass, force, pressure, density, allocMemoryType);
	}

	///Destructor.
	virtual ~Particles(){}

	virtual unsigned int size() const {return m_pos->size();}

	///Create children.
	void createChildren(bool ppPos, bool pPos, bool velocity, bool acceleration, bool mass, bool force, bool pressure, bool density,
		unsigned int size=0, MemoryType allocMemoryType=BufferSet::HOST);

public:

	///Particle prev prev position.
	Points* m_ppPos;

	///Particle prev position.
	Points* m_pPos;

	///Particle position. Always exists.
	Points* m_pos;

	///Particle velocity.
	Points* m_velocity;

	///Particle acceleration.
	Points* m_acceleration;

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
