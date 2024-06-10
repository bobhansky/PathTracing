#pragma once

#include "IIntersectStrategy.h"
#include "BVH.hpp"
#include "Vector.hpp"

class BVHStrategy : public IIntersectStrategy {
	virtual void UpdateInter(Intersection& inter, Scene& sce, 
		const Vector3d & rayOrig, const Vector3d& rayDir)override {
		inter = getIntersection(sce.BVHaccelerator->getNode(), rayOrig, rayDir);
	}

	virtual double getShadowCoeffi(Scene& sce, Intersection& p, Vector3d& lightPos) override{
		Vector3d orig = p.pos;
		orig = orig + 0.0005f * p.nDir;
		Vector3d raydir = normalized(lightPos - orig);

		double distance = (lightPos - orig).norm();
		
		return ShadowHelper(sce.BVHaccelerator->getNode(), orig, raydir, distance);
	}

	double ShadowHelper(BVHNode* node, const Vector3d& rayOrig, const Vector3d& rayDir, double dis) {
		if (!node) return 1;
		// if ray miss this bound
		if (!node->bound.IntersectRay(rayOrig, rayDir))
			return 1;

		Intersection inter;
		// if the node is a leaf node
		// then test the intersection of ray and objects
		// in this project, leaf box always contain only 1 object
		if (!node->left && !node->right) {
			node->obj->intersect(rayOrig, rayDir, inter);
			if (inter.intersected && inter.t < dis)
				return (1-inter.mtlcolor.alpha);
			return 1;
		}

		// if node is a internal node
		double l = ShadowHelper(node->left, rayOrig, rayDir, dis);
		double r = ShadowHelper(node->right, rayOrig, rayDir, dis);

		return l * r;
	}


};