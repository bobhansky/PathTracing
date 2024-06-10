#pragma once

#include <float.h>
#include "global.hpp"
#include "Vector.hpp"
#include "Material.hpp"
//#include "Object.hpp"


class Object;	// circular dependency issue
				// see https://stackoverflow.com/questions/23283080/compiler-error-c4430-missing-type-specifier-int-assumed

class Intersection {
public:

	bool intersected = false;
	double t = FLT_MAX;	// pos = rayPos + t * rayDir
	Vector3d pos;
	Vector3d nDir; // normal direction

	Vector2d textPos;	// texture coordinates if any	(-1, -1) means no texture
	int diffuseIndex = -1;
	int normalMapIndex = -1;
	int roughnessMapIndex = -1;
	int metallicMapIndex = -1;


	Material mtlcolor;
	Object *obj = nullptr;		// this intersection is on which object	

};