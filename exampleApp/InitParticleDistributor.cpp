#include "InitParticleDistributor.h"

using namespace hohehohe2;


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void InitParticleDistributor::set(PointSet& pos, PointSet& velocity, float restLength, unsigned int id)
{
	switch (id)
	{
	case 0:
		placement0_(pos, velocity, restLength);
		break;
	case 1:
		placement1_(pos, velocity, restLength);
		break;
	case 2:
		placement2_(pos, velocity, restLength);
		break;
	case 3:
		placement3_(pos, velocity, restLength);
		break;
	case 4:
		placement4_(pos, velocity, restLength);
		break;
	default:
		placement0_(pos, velocity, restLength);
		break;
	}

	pos.setClean(HOST);
	velocity.setClean(HOST);
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void InitParticleDistributor::placement0_(PointSet& pos, PointSet& velocity, float restLength)
{
	pos.setSize(1);
	pos.allocate(HOST);
	velocity.setSize(1);
	velocity.allocate(HOST);

	float* pxs = pos.xs(true);
	float* pys = pos.ys(true);
	float* pzs = pos.zs(true);
	float* vxs = velocity.xs(true);
	float* vys = velocity.ys(true);
	float* vzs = velocity.zs(true);

	pxs[0] = 0.0f;
	pys[0] = 0.0f;
	pzs[0] = 0.0f;
	vxs[0] = 0.0f;
	vys[0] = 0.0f;
	vzs[0] = 0.0f;
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void InitParticleDistributor::placement1_(PointSet& pos, PointSet& velocity, float restLength)
{
	pos.setSize(2);
	pos.allocate(HOST);
	velocity.setSize(2);
	velocity.allocate(HOST);

	float* pxs = pos.xs(true);
	float* pys = pos.ys(true);
	float* pzs = pos.zs(true);
	float* vxs = velocity.xs(true);
	float* vys = velocity.ys(true);
	float* vzs = velocity.zs(true);

	pxs[0] = -restLength;
	pys[0] = 0.0f;
	pzs[0] = 0.0f;
	vxs[0] = 0.0f;
	vys[0] = 0.0f;
	vzs[0] = 0.0f;

	pxs[1] = restLength;
	pys[1] = 0.0f;
	pzs[1] = 0.0f;
	vxs[1] = 0.0f;
	vys[1] = 0.0f;
	vzs[1] = 0.0f;
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void InitParticleDistributor::placement2_(PointSet& pos, PointSet& velocity, float restLength)
{
	const unsigned int NUM_LINES = 10;
	pos.setSize(NUM_LINES * NUM_LINES * NUM_LINES);
	pos.allocate(HOST);
	velocity.setSize(NUM_LINES * NUM_LINES * NUM_LINES);
	velocity.allocate(HOST);

	float* pxs = pos.xs(true);
	float* pys = pos.ys(true);
	float* pzs = pos.zs(true);
	float* vxs = velocity.xs(true);
	float* vys = velocity.ys(true);
	float* vzs = velocity.zs(true);

	for (unsigned int i = 0; i < NUM_LINES; ++i)
	{
		for (unsigned int j = 0; j < NUM_LINES; ++j)
		{
			for (unsigned int k = 0; k < NUM_LINES; ++k)
			{
				unsigned int pid = i * NUM_LINES * NUM_LINES + j * NUM_LINES + k;
				pxs[pid] = restLength / 1.1f * (i - NUM_LINES / 2);
				pys[pid] = restLength / 1.1f * (j - NUM_LINES / 2);
				pzs[pid] = restLength / 1.1f * (k - NUM_LINES / 2);;
				vxs[pid] = 0.0f;
				vys[pid] = 0.0f;
				vzs[pid] = 0.0f;
			}
		}
	}
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void InitParticleDistributor::placement3_(PointSet& pos, PointSet& velocity, float restLength)
{
	pos.setSize(10);
	pos.allocate(HOST);
	velocity.setSize(10);
	velocity.allocate(HOST);

	float* pxs = pos.xs(true);
	float* pys = pos.ys(true);
	float* pzs = pos.zs(true);
	float* vxs = velocity.xs(true);
	float* vys = velocity.ys(true);
	float* vzs = velocity.zs(true);

	for (unsigned int i = 0; i < 10; ++i)
	{
		pxs[i] = 0;
		pys[i] = restLength * (9 - i);
		pzs[i] = 0;
		vxs[i] = 0.0f;
		vys[i] = 0.0f;
		vzs[i] = 0.0f;
	}
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void InitParticleDistributor::placement4_(PointSet& pos, PointSet& velocity, float restLength)
{
	const unsigned int NUM_LINES = 15;
	pos.setSize(NUM_LINES * NUM_LINES * NUM_LINES);
	pos.allocate(HOST);
	velocity.setSize(NUM_LINES * NUM_LINES * NUM_LINES);
	velocity.allocate(HOST);

	float* pxs = pos.xs(true);
	float* pys = pos.ys(true);
	float* pzs = pos.zs(true);
	float* vxs = velocity.xs(true);
	float* vys = velocity.ys(true);
	float* vzs = velocity.zs(true);

	for (unsigned int i = 0; i < NUM_LINES; ++i)
	{
		for (unsigned int j = 0; j < NUM_LINES; ++j)
		{
			for (unsigned int k = 0; k < NUM_LINES; ++k)
			{
				unsigned int pid = i * NUM_LINES * NUM_LINES + j * NUM_LINES + k;
				pxs[pid] = restLength / 1.0f * (i - NUM_LINES / 2);
				pys[pid] = restLength / 1.0f * (j - NUM_LINES / 2);
				pzs[pid] = restLength / 1.0f * (k - NUM_LINES / 2);;
				vxs[pid] = 0.0f;
				vys[pid] = 0.0f;
				vzs[pid] = 0.0f;
			}
		}
	}
}
