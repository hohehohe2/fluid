#ifndef hohe_Particles_H
#define hohe_Particles_H

#include <cudaCommon/Buffer.h>
#include "basicGeos.h"

namespace hohehohe2
{

//-------------------------------------------------------------------
//-------------------------------------------------------------------
///Particle cloud.
class Particles : public BufferSetSameSizedMember
{

public:

	///Constructor.
	/**
	@param hasVelocity Set true if the partcles have a velocity buffer.
	@param hasMass Set true if the partcles have a mass buffer.
	@param hasForce Set true if the partcles have a force buffer.
	@param size Number of particles.
	**/
	Particles(bool hasVelocity, bool hasMass, bool hasForce, unsigned int size=0);


	///Destructor.
	~Particles(){}

	virtual unsigned int size() const {return m_pos.size();}

public:

	///Particle prev prev position.
	Points m_ppPos;

	///Particle prev position.
	Points m_pPos;

	///Particle position.
	Points m_pos;

	///Particle velocity.
	Buffer < float > m_velocity;

	///Particle mass.
	Buffer < float > m_mass;

	///Particle force.
	Points m_force;

};

}

#endif
