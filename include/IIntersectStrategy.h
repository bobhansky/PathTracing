#pragma once

#include "Scene.hpp"
#include "Intersection.hpp"

// intersection test can be the basic way or the BVH way
class IIntersectStrategy {
public:
	// intersection with object
	virtual void UpdateInter(Intersection& inter, Scene& sce,
		const Vector3d& rayOrig, const Vector3d& rayDir) = 0;

	// calculate shadow coefficient, for hard shadow
	virtual double getShadowCoeffi(Scene& sce, Intersection& p, Vector3d& lightpos) = 0;


};