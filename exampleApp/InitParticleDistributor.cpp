#include "InitParticleDistributor.h"
#include <sph/ParticlesFluid.h>
#include <sph/ParticlesWall.h>

using namespace hohehohe2;


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void InitParticleDistributor::set(PointSet& pos, PointSet& velocity, PointSet& posWall, float equilibriumDistance, unsigned int id)
{
	switch (id)
	{
	case 0:
		placement0_(pos, velocity, posWall, equilibriumDistance);
		break;
	case 1:
		placement1_(pos, velocity, posWall, equilibriumDistance);
		break;
	case 2:
		placement2_(pos, velocity, posWall, equilibriumDistance);
		break;
	case 3:
		placement3_(pos, velocity, posWall, equilibriumDistance);
		break;
	case 4:
		placement4_(pos, velocity, posWall, equilibriumDistance);
		break;
	case 5:
		placement5_(pos, velocity, posWall, equilibriumDistance);
		break;
	case 6:
		placement6_(pos, velocity, posWall, equilibriumDistance);
		break;
	case 7:
		placement7_(pos, velocity, posWall, equilibriumDistance);
		break;
	default:
		placement0_(pos, velocity, posWall, equilibriumDistance);
		break;
	}

	pos.setClean(HOST);
	velocity.setClean(HOST);
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void InitParticleDistributor::placement0_(PointSet& pos, PointSet& velocity, PointSet& posWall, float equilibriumDistance)
{
	pos.setSize(1);
	pos.allocate(HOST);
	velocity.setSize(1);
	velocity.allocate(HOST);

	float* pxs = pos.xs(HOST);
	float* pys = pos.ys(HOST);
	float* pzs = pos.zs(HOST);
	float* vxs = velocity.xs(HOST);
	float* vys = velocity.ys(HOST);
	float* vzs = velocity.zs(HOST);

	pxs[0] = 0.0f;
	pys[0] = 0.0f;
	pzs[0] = 0.0f;
	vxs[0] = 0.0f;
	vys[0] = 0.0f;
	vzs[0] = 0.0f;

	pos.setClean(HOST);
	velocity.setClean(HOST);
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void InitParticleDistributor::placement1_(PointSet& pos, PointSet& velocity, PointSet& posWall, float equilibriumDistance)
{
	pos.setSize(2);
	pos.allocate(HOST);
	velocity.setSize(2);
	velocity.allocate(HOST);

	float* pxs = pos.xs(HOST);
	float* pys = pos.ys(HOST);
	float* pzs = pos.zs(HOST);
	float* vxs = velocity.xs(HOST);
	float* vys = velocity.ys(HOST);
	float* vzs = velocity.zs(HOST);

	pxs[0] = -equilibriumDistance / 2;
	pys[0] = 0.0f;
	pzs[0] = 0.0f;
	vxs[0] = 0.0f;
	vys[0] = 0.0f;
	vzs[0] = 0.0f;

	pxs[1] = equilibriumDistance / 2;
	pys[1] = 0.0f;
	pzs[1] = 0.0f;
	vxs[1] = 0.0f;
	vys[1] = 0.0f;
	vzs[1] = 0.0f;

	pos.setClean(HOST);
	velocity.setClean(HOST);
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void InitParticleDistributor::placement2_(PointSet& pos, PointSet& velocity, PointSet& posWall, float equilibriumDistance)
{
	const int NUM_LINES = 10;
	pos.setSize(NUM_LINES * NUM_LINES * NUM_LINES);
	pos.allocate(HOST);
	velocity.setSize(NUM_LINES * NUM_LINES * NUM_LINES);
	velocity.allocate(HOST);

	float* pxs = pos.xs(HOST);
	float* pys = pos.ys(HOST);
	float* pzs = pos.zs(HOST);
	float* vxs = velocity.xs(HOST);
	float* vys = velocity.ys(HOST);
	float* vzs = velocity.zs(HOST);

	for (int i = 0; i < NUM_LINES; ++i)
	{
		for (int j = 0; j < NUM_LINES; ++j)
		{
			for (int k = 0; k < NUM_LINES; ++k)
			{
				unsigned int pid = i * NUM_LINES * NUM_LINES + j * NUM_LINES + k;
				pxs[pid] = equilibriumDistance / 1.0f * (i - NUM_LINES / 2);
				pys[pid] = equilibriumDistance / 1.0f * j;
				pzs[pid] = equilibriumDistance / 1.0f * (k - NUM_LINES / 2);
				vxs[pid] = 0.0f;
				vys[pid] = 0.0f;
				vzs[pid] = 0.0f;
			}
		}
	}

	pos.setClean(HOST);
	velocity.setClean(HOST);
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void InitParticleDistributor::placement3_(PointSet& pos, PointSet& velocity, PointSet& posWall, float equilibriumDistance)
{
	pos.setSize(10);
	pos.allocate(HOST);
	velocity.setSize(10);
	velocity.allocate(HOST);

	float* pxs = pos.xs(HOST);
	float* pys = pos.ys(HOST);
	float* pzs = pos.zs(HOST);
	float* vxs = velocity.xs(HOST);
	float* vys = velocity.ys(HOST);
	float* vzs = velocity.zs(HOST);

	for(int i = 0; i < 10; ++i)
	{
		pxs[i] = 0;
		pys[i] = equilibriumDistance * (9 - i);
		pzs[i] = 0;
		vxs[i] = 0.0f;
		vys[i] = 0.0f;
		vzs[i] = 0.0f;
	}

	pos.setClean(HOST);
	velocity.setClean(HOST);
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void InitParticleDistributor::placement4_(PointSet& pos, PointSet& velocity, PointSet& posWall, float equilibriumDistance)
{
	const int NUM_LINES = 50;
	pos.setSize(NUM_LINES * NUM_LINES * NUM_LINES);
	pos.allocate(HOST);
	velocity.setSize(NUM_LINES * NUM_LINES * NUM_LINES);
	velocity.allocate(HOST);

	float* pxs = pos.xs(HOST);
	float* pys = pos.ys(HOST);
	float* pzs = pos.zs(HOST);
	float* vxs = velocity.xs(HOST);
	float* vys = velocity.ys(HOST);
	float* vzs = velocity.zs(HOST);

	for(int i = 0; i < NUM_LINES; ++i)
	{
		for(int j = 0; j < NUM_LINES; ++j)
		{
			for(int k = 0; k < NUM_LINES; ++k)
			{
				unsigned int pid = i * NUM_LINES * NUM_LINES + j * NUM_LINES + k;
				pxs[pid] = equilibriumDistance / 1.0f * (i - NUM_LINES / 2);
				pys[pid] = equilibriumDistance / 1.0f * (j - NUM_LINES / 2);
				pzs[pid] = equilibriumDistance / 1.0f * (k - NUM_LINES / 2);;
				vxs[pid] = 0.0f;
				vys[pid] = 0.0f;
				vzs[pid] = 0.0f;
			}
		}
	}

	pos.setClean(HOST);
	velocity.setClean(HOST);
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void InitParticleDistributor::placement5_(PointSet& pos, PointSet& velocity, PointSet& posWall, float equilibriumDistance)
{
	const int NUM_LINES = 10;
	pos.setSize(NUM_LINES * NUM_LINES * NUM_LINES);
	pos.allocate(HOST);
	velocity.setSize(NUM_LINES * NUM_LINES * NUM_LINES);
	velocity.allocate(HOST);

	float* pxs = pos.xs(HOST);
	float* pys = pos.ys(HOST);
	float* pzs = pos.zs(HOST);
	float* vxs = velocity.xs(HOST);
	float* vys = velocity.ys(HOST);
	float* vzs = velocity.zs(HOST);

	for (int i = 0; i < NUM_LINES; ++i)
	{
		for (int j = 0; j < NUM_LINES; ++j)
		{
			for (int k = 0; k < NUM_LINES; ++k)
			{
				unsigned int pid = i * NUM_LINES * NUM_LINES + j * NUM_LINES + k;
				pxs[pid] = equilibriumDistance / 1.1f * (i - NUM_LINES / 2);
				pys[pid] = equilibriumDistance / 1.1f * j;
				pzs[pid] = equilibriumDistance / 1.1f * (k - NUM_LINES / 2);
				vxs[pid] = 0.0f;
				vys[pid] = 0.0f;
				vzs[pid] = 0.0f;
			}
		}
	}

	pos.setClean(HOST);
	velocity.setClean(HOST);

	const int NUM_GROUND_LINES = 100;
	posWall.setSize(NUM_GROUND_LINES * NUM_GROUND_LINES);
	posWall.allocate(HOST);

	float* wpxs = posWall.xs(HOST);
	float* wpys = posWall.ys(HOST);
	float* wpzs = posWall.zs(HOST);

	for(int i = 0; i < NUM_GROUND_LINES; ++i)
	{
		for(int k = 0; k < NUM_GROUND_LINES; ++k)
		{
			unsigned int pid = i * NUM_GROUND_LINES + k;
			wpxs[pid] = equilibriumDistance * (i - NUM_GROUND_LINES / 2);
			wpys[pid] = -equilibriumDistance * 10;
			wpzs[pid] = equilibriumDistance * (k - NUM_GROUND_LINES / 2);
		}
	}

	posWall.setClean(HOST);
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void InitParticleDistributor::placement6_(PointSet& pos, PointSet& velocity, PointSet& posWall, float equilibriumDistance)
{
	pos.setSize(4);
	pos.allocate(HOST);
	velocity.setSize(4);
	velocity.allocate(HOST);

	float* pxs = pos.xs(HOST);
	float* pys = pos.ys(HOST);
	float* pzs = pos.zs(HOST);
	float* vxs = velocity.xs(HOST);
	float* vys = velocity.ys(HOST);
	float* vzs = velocity.zs(HOST);

	pxs[0] = 0.0;
	pys[0] = 0.0f;
	pzs[0] = 0.0f;
	vxs[0] = 0.0f;
	vys[0] = 0.0f;
	vzs[0] = 0.0f;

	pxs[1] = equilibriumDistance;
	pys[1] = 0.0f;
	pzs[1] = 0.0f;
	vxs[1] = 0.0f;
	vys[1] = 0.0f;
	vzs[1] = 0.0f;

	pxs[2] = 0.0f;
	pys[2] = equilibriumDistance;
	pzs[2] = 0.0f;
	vxs[2] = 0.0f;
	vys[2] = 0.0f;
	vzs[2] = 0.0f;

	pxs[3] = equilibriumDistance;
	pys[3] = equilibriumDistance;
	pzs[3] = 0.0f;
	vxs[3] = 0.0f;
	vys[3] = 0.0f;
	vzs[3] = 0.0f;

	pos.setClean(HOST);
	velocity.setClean(HOST);
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
void InitParticleDistributor::placement7_(PointSet& pos, PointSet& velocity, PointSet& posWall, float equilibriumDistance)
{
	pos.setSize(1);
	pos.allocate(HOST);
	velocity.setSize(1);
	velocity.allocate(HOST);

	float* pxs = pos.xs(HOST);
	float* pys = pos.ys(HOST);
	float* pzs = pos.zs(HOST);
	float* vxs = velocity.xs(HOST);
	float* vys = velocity.ys(HOST);
	float* vzs = velocity.zs(HOST);

	pxs[0] = 0.0f;
	pys[0] = 0.0f;
	pzs[0] = 0.0f;
	vxs[0] = 0.0f;
	vys[0] = 0.0f;
	vzs[0] = 0.0f;

	pos.setClean(HOST);
	velocity.setClean(HOST);

	const int NUM_GROUND_LINES = 100;
	posWall.setSize(NUM_GROUND_LINES * NUM_GROUND_LINES);
	posWall.allocate(HOST);

	float* wpxs = posWall.xs(HOST);
	float* wpys = posWall.ys(HOST);
	float* wpzs = posWall.zs(HOST);

	for(int i = 0; i < NUM_GROUND_LINES; ++i)
	{
		for(int k = 0; k < NUM_GROUND_LINES; ++k)
		{
			unsigned int pid = i * NUM_GROUND_LINES + k;
			wpxs[pid] = equilibriumDistance * (i - NUM_GROUND_LINES / 2);
			wpys[pid] = -equilibriumDistance * 10;
			wpzs[pid] = equilibriumDistance * (k - NUM_GROUND_LINES / 2);
		}
	}

	posWall.setClean(HOST);
}