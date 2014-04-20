#include "Particles.h"

using namespace hohehohe2;


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void Particles::createChildren(bool ppPos, bool pPos, bool velocity, bool acceleration, bool mass, bool force, bool pressure, bool density,
	unsigned int size, MemoryType allocMemoryType)
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
		addChild(m_velocity = new Points("particle velocity"));
	}
	if (acceleration)
	{
		removeChild(m_acceleration);
		addChild(m_acceleration = new Points("particle acceleration"));
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
