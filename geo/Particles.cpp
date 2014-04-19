#include "Particles.h"

using namespace hohehohe2;

//-------------------------------------------------------------------
//-------------------------------------------------------------------
void Particles::init(unsigned int size, bool pppos, bool ppos, bool pos, bool velocity, bool mass, bool force, bool pressure, bool density)
{
	clearChild();
	if (pppos)
	{
		m_ppPos.reset(new Points("particles pppos"));
		addChild(m_ppPos);
	}
	if (ppos)
	{
		m_pPos.reset(new Points("particles ppos"));
		addChild(m_pPos);
	}
	if (pos)
	{
		m_pos.reset(new Points("particles pos"));
		addChild(m_pos);
	}
	if (velocity)
	{
		m_velocity.reset(new BufferFloat("particles velocity"));
		addChild(m_velocity);
	}
	if (mass)
	{
		m_mass.reset(new BufferFloat("particles mass"));
		addChild(m_mass);
	}
	if (force)
	{
		m_force.reset(new Points("particles force"));
		addChild(m_force);
	}
	if (pressure)
	{
		m_pressure.reset(new BufferFloat("particles pressure"));
		addChild(m_pressure);
	}
	if (density)
	{
		m_density.reset(new BufferFloat("particles density"));
		addChild(m_density);
	}
	setSize(size);
}
