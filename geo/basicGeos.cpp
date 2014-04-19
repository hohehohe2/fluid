#include "basicGeos.h"

using namespace hohehohe2;

//-------------------------------------------------------------------
//-------------------------------------------------------------------
void Points::init(unsigned int size, bool xs, bool ys, bool zs)
{
	clearChild();
	if (xs)
	{
		m_xs.reset(new BufferFloat("point xs"));
		addChild(m_xs);
	}
	if (ys)
	{
		m_ys.reset(new BufferFloat("point ys"));
		addChild(m_ys);
	}
	if (zs)
	{
		m_ys.reset(new BufferFloat("point zs"));
		addChild(m_zs);
	}
	setSize(size);
}



//-------------------------------------------------------------------
//-------------------------------------------------------------------
void Lines::init(unsigned int size, bool starts, bool ends)
{
	clearChild();
	if (starts)
	{
		m_starts.reset(new Points("line starts"));
		addChild(m_starts);
	}
	if (ends)
	{
		m_ends.reset(new Points("line ends"));
		addChild(m_ends);
	}
	setSize(size);
}



//-------------------------------------------------------------------
//-------------------------------------------------------------------
void Triangles::init(unsigned int size, bool pos, bool v0s, bool v1s, bool v2s)
{
	clearChild();
	if (pos)
	{
		m_pos.reset(new Points("triangle pos"));
		addChild(m_pos);
	}
	if (v0s)
	{
		m_v0s.reset(new BufferUInt("triangle v0s"));
		addChild(m_v0s);
	}
	if (v1s)
	{
		m_v1s.reset(new BufferUInt("triangle v1s"));
		addChild(m_v1s);
	}
	if (v2s)
	{
		m_v2s.reset(new BufferUInt("triangle v2s"));
		addChild(m_v2s);
	}
	setSize(size);
}
