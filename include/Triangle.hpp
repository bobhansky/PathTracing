#pragma once

#include "Object.hpp"
#include "Intersection.hpp"



class Triangle : public Object{
public:
	// 3 vertices, from lower left counterclockwise 
	Vector3d v0;
	Vector3d v1;
	Vector3d v2;

	// normal direction of 3 vertices
	Vector3d n0, n1, n2;
	Vector2d uv0, uv1, uv2;	// texture coordinate u,v   (-1,-1) initially means no texture 

	// Using Moller Trumbore Algorithm to update the intersection between ray and triangle
	// solve with Cramer's rule
	// https://www.geeksforgeeks.org/system-linear-equations-three-variables-using-cramers-rule/#
	// i didnt read it carefully
	bool intersect(const Vector3d& orig, const Vector3d& dir, Intersection& inter) override {
		
		Vector3d E1 = v1 - v0;
		Vector3d E2 = v2 - v0;          // v2 - v1   get the strange res
		Vector3d S = orig - v0;
		Vector3d S1 = crossProduct(dir, E2);		// pvec
		Vector3d S2 = crossProduct(S, E1);			

		// check if ray parallel to the surface
		// and 
		// under surface ray					
		Vector3d normal = crossProduct(E1, E2);
		normal = normalized(normal);

		// if the ray is parallel to the surface, then no inter
		// I CONSIDER undersurface ray, so it's not >=
		if (FLOAT_EQUAL( dir.dot(normal), 0.0))
			return false;
		

		Vector3d rightVec(S2.dot(E2), S1.dot(S),S2.dot(dir));	// right handside part in my note games101
		if (S1.dot(E1) == 0.0) return false;
		double left = 1.0f / S1.dot(E1);

		Vector3d res = left * rightVec;		// t u v are in res now

		if (res.x  > 0 && 1 - res.y - res.z  > 0 && res.y  > 0 && res.z  > 0) {
			// then update intersection
			inter.intersected = true;
			inter.obj = this;
			inter.t = res.x;
			inter.pos = orig + inter.t * dir;
			inter.mtlcolor = this->mtlcolor;
			inter.nDir = normalized((n0 * (1 - res.y - res.z)) + n1 * res.y + n2 * res.z);	// smooth shading

			// inter.nDir = normalized(crossProduct(E1, E2));	// flat shading

			// calculate texture coordinates
			if (isTextureActivated)
			{
				inter.textPos = uv0 * (1 - res.y - res.z) + uv1 * res.y + uv2 * res.z;
				inter.diffuseIndex = this->textureIndex;
				inter.normalMapIndex = this->normalMapIndex;
				inter.roughnessMapIndex = roughnessMapIndex;
				inter.metallicMapIndex = metallicMapIndex;
			}
			return true;
		}
		return false;
		

	
	}


	// get the alpha beta gamma in barycentric corrdinate
	// update a b g in the parameter list
	// return true if we can find a b g
	// false otherwise
	bool getBarycentric(const Vector3d& point,
		double& alpha, double& beta, double& gamma) {
		// triangles.pdf   page 69

		Vector3d e1 = v1 - v0;
		Vector3d e2 = v2 - v0;
		Vector3d ep = point - v0;

		double d11 = e1.dot(e1);
		double d12 = e1.dot(e2);
		double d22 = e2.dot(e2);
		double d1p = e1.dot(ep);
		double d2p = e2.dot(ep);

		double det = d11 * d22 - d12 * d12;
		if (FLOAT_EQUAL(0.0, det)) return false;		// 3 vertices on a line: not a triangle

		beta = (d22 * d1p - d12 * d2p) / det;
		gamma = (d11 * d2p - d12 * d1p) / det;
		alpha = 1 - beta - gamma;
		return true;
	}

	void initializeBound() override {
		bound = BoundBox(v0, v1);
		bound = Union(bound, v2);
	}

	double getArea() override {
		Vector3d e1 = v1 - v0;
		Vector3d e2 = v2 - v0;

		e1 = v1 - v0;
		e2 = v2 - v0;
		return  crossProduct(e1, e2).norm() * 0.5f;
	}

	void samplePoint(Intersection& inter, double& pdf) override {
		
		double u = getRandomFloat();
		double v = getRandomFloat() * (1 - u);

		Vector3d pos = (1 - u - v) * v0 + u * v1 + v * v2;


		inter.pos = pos;
		inter.nDir = (1 - u - v) * n0 + u * n1 + v * n2;
		inter.nDir = normalized(inter.nDir);
		inter.intersected = true;
		inter.mtlcolor = mtlcolor;
		inter.obj = this;
		

		double area = getArea();
		pdf = 1.0 / area ;

		
		if (isTextureActivated) {
			inter.normalMapIndex = normalMapIndex;
			inter.textPos = uv0 * (1 - u - v) + uv1 * u + uv2 * v;
			inter.diffuseIndex = textureIndex;
		}
	}
};