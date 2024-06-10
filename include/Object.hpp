#pragma once

#include "Vector.hpp"
#include "global.hpp"
#include "Material.hpp"
#include "Intersection.hpp"
#include "BoundBox.hpp"

enum OBJTYPE
{
	TRIANGLE,
	SPEHRE
};

class Object {
public:
	virtual ~Object() {};
	// check if ray intersect with this object
	// if intersect then update inter data
	// orig: ray origin
	// dir: ray direction
	virtual bool intersect(const Vector3d& orig, const Vector3d& dir, Intersection& inter) = 0;

	Object() {

	};

	OBJTYPE objectType;
	Material mtlcolor;

	bool isTextureActivated = false;
	int textureIndex = -1;		// diffuse color map index
	int normalMapIndex = -1;    // normal map index
	int roughnessMapIndex = -1;
	int metallicMapIndex = -1;

	BoundBox bound;
	// initialize the bound of this object
	virtual void initializeBound() = 0;
	// get the surface area of this object
	virtual double getArea() = 0;
	// randomly sample a point on the surface of this object
	virtual void samplePoint(Intersection& inter, double& pdf) = 0;
	
};