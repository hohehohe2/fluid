#include "Particles.h"

using namespace hohehohe2;

//-------------------------------------------------------------------
//-------------------------------------------------------------------
Particles::Particles(bool hasVelocity, bool hasMass, bool hasForce, unsigned int size) : m_pos(size)
{
	addMember_(m_pos);
	if (hasVelocity)
	{
		addMember_(m_velocity);
		m_velocity.setSize(size);
	}
	if (hasMass)
	{
		addMember_(m_mass);
		m_mass.setSize(size);
	}
	if (hasForce)
	{
		addMember_(m_force);
		m_force.setSize(size);
	}
}
