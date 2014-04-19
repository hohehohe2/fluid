#ifndef hohe_Particles_H
#define hohe_Particles_H

#include <cudaCommon/Buffer.h>
#include "basicGeos.h"

namespace hohehohe2
{

//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Particle cloud.
class Particles : public BufferSetSameSized
{

public:

	typedef std::shared_ptr < Particles > SPtr;

	///Constructor.
	Particles(const std::string& name="") : BufferSetSameSized(name){}

	///Destructor.
	virtual ~Particles(){}

	virtual unsigned int size() const {return m_pos->size();}

	void init(unsigned int size, bool pppos, bool ppos, bool pos, bool velocity, bool mass, bool force, bool pressure, bool density);

public:

	///Particle prev prev position.
	Points::SPtr m_ppPos;

	///Particle prev position.
	Points::SPtr m_pPos;

	///Particle position. Always exists.
	Points::SPtr m_pos;

	///Particle velocity.
	BufferFloat::SPtr m_velocity;

	///Particle mass.
	BufferFloat::SPtr m_mass;

	///Particle force.
	Points::SPtr m_force;

	///Particle pressure.
	BufferFloat::SPtr m_pressure;

	///Particle density.
	BufferFloat::SPtr m_density;

};

}

#endif
