#pragma once

#include "Object.hpp"


class Sphere : public Object{
public:
	Vector3d centerPos;		// center position
	double radius;

	Sphere() {
		centerPos = Vector3d(0.0, 0.0, 0.0);
		radius = 1.0;
	}

	Sphere(double x, double y, double z, double r) {
		centerPos = Vector3d(x, y, z);
		radius = r;
	}


	virtual ~Sphere() {};

	// check if the ray will intersect with this sphere or not
	// if true, then set the nearest time of intersection as tNear, and mtlColor
	bool intersect(const Vector3d& orig, const Vector3d& dir, Intersection& inter) override {
		double A = 1.0;	// here is 1 since we are using normalized vector
		double B = 2 * (dir.x * (orig.x - centerPos.x) + dir.y * (orig.y - centerPos.y)
			+ dir.z * (orig.z - centerPos.z));
		double C = pow((orig.x - centerPos.x), 2) + pow((orig.y - centerPos.y), 2)
			+ pow((orig.z - centerPos.z), 2) - radius * radius;

		double t1 = 0;
		double t2 = 0;
		solveQuadratic(t1, t2, A, B, C);

		inter.intersected = false;

		// miss, no real solution
		if (FLOAT_EQUAL(t1, FLT_MAX) && FLOAT_EQUAL(t2, FLT_MAX)) {
			return false;
		}

		// ONE soluion
		else if (FLOAT_EQUAL(t1, t2)) {
			// intersection is behind the ray direction, then false
			if (t1 < 0 ) return false;
			else {
				inter.t = t1;
				// update intersection data
				inter.intersected = true;
				inter.obj = this;
				inter.mtlcolor = this->mtlcolor;
				inter.pos = orig + inter.t * dir;
				inter.nDir = normalized(inter.pos - centerPos);
				if (isTextureActivated)
				{
					// calculate the uv coordinate of this intersection
					double u, v;
					double phi = acos(inter.nDir.z);
					v = phi / M_PI;

					double theta = atan2(inter.nDir.y, inter.nDir.x);	// return [-pi, pi]
					// we need to map it to [0, 1]
					if (theta < 0) theta += 2 * M_PI;	// trigonometric functions are periodic
					u = (theta / (2.0 * M_PI));			// 0 + [0, 1]    then if theta == 0, it is the left most point in width

					// or 
					// u = 0.5 + (theta / (2.0 * M_PI));  // 0.5 + [-0.5, 0.5]	  then if theta == 0, it is the middle point in width

					inter.textPos = Vector2d(u, v);
					inter.diffuseIndex = this->textureIndex;
					inter.normalMapIndex = normalMapIndex;
					inter.roughnessMapIndex = roughnessMapIndex;
					inter.metallicMapIndex = metallicMapIndex;
				}
				return true;
			}

		}
		// TWO solution
		else {
			if (t1 > 0 && t2 > 0) {
				inter.t = t1;
			}
			else if (t1 > 0 && t2 < 0) {
				inter.t = t1;
			}
			else if (t1 < 0 && t2 > 0) {
				inter.t = t2;
			}
			else return false;

			// update intersection data
			inter.intersected = true;
			inter.obj = this;
			inter.mtlcolor = this->mtlcolor;
			inter.pos = orig + inter.t * dir;
			inter.nDir = normalized(inter.pos - centerPos);
			if (isTextureActivated)
			{
				// calculate the uv coordinate of this intersection
				double u, v;
				double phi = acos(inter.nDir.z);	// return [0, pi]
				v = phi / M_PI;

				double theta = atan2(inter.nDir.y, inter.nDir.x);	// return [-pi, pi]
				// we need to map it to [0, 1]
				if (theta < 0) theta += 2 * M_PI;	// trigonometric functions are periodic
				u = (theta / (2.0 * M_PI));			// 0 + [0, 1]    then if theta == 0, it is the left most point in width

				// or 
				// u = 0.5 + (theta / (2.0 * M_PI));  // 0.5 + [-0.5, 0.5]	  then if theta == 0, it is the middle point in width

				inter.textPos = Vector2d(u, v);
				inter.diffuseIndex = this->textureIndex;
				inter.normalMapIndex = normalMapIndex;
				inter.roughnessMapIndex = roughnessMapIndex;
				inter.metallicMapIndex = metallicMapIndex;
			}
			return true;
		}
		return false;
	}


	void initializeBound() override {
		Vector3d min = { centerPos.x - radius, centerPos.y - radius, centerPos.z - radius };
		Vector3d max = { centerPos.x + radius, centerPos.y + radius, centerPos.z + radius };
		bound = BoundBox(min, max);
	}

	double getArea() override {
		return radius * radius * M_PI;
	}

	void samplePoint(Intersection& inter, double& pdf) override {
		// 10/1/2023:
		// is it correct? maybe I should Inverse Trasnform sampling?

		// theta: angle from +x to yz plane	
		// phi: angle from +z to xy plane
		// 7/6/2023: 

		double theta = getRandomFloat() * 2 * M_PI;
		double phi = getRandomFloat() * M_PI;

		inter.pos.x = centerPos.x + radius * cos(theta) * sin(phi);
		inter.pos.y = centerPos.y + radius * sin(theta) * sin(phi);
		inter.pos.z = centerPos.z + radius * cos(phi);

		inter.nDir = normalized(inter.pos - centerPos);

		inter.intersected = true;
		inter.mtlcolor = mtlcolor;
		inter.obj = this;

		// I don't update the texture corrdinate here	7/2/2023

		pdf = 1 / getArea();
	}
};