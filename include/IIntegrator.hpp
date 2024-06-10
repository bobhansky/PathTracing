#pragma once

#include<mutex>

#include "Vector.hpp"
#include "global.hpp"
#include "PPMGenerator.hpp"
#include "IIntersectStrategy.h"

std::mutex sampleLight_mutex;
std::vector<std::string> records;	// ray information, thread independent string

class IIntegrator {
public:
	virtual void integrate(PPMGenerator* g) = 0;

protected:
	IIntersectStrategy* interStrategy;
	PPMGenerator* g;
};



void changeNormalDir(Intersection& inter, PPMGenerator* g) {
	Texture* nMap = g->normalMaps.at(inter.normalMapIndex);
	Vector3d color = nMap->getRGBat(inter.textPos.x, inter.textPos.y);

	switch (inter.obj->objectType)
	{
	case TRIANGLE: {
		Triangle* t = static_cast<Triangle*>(inter.obj);
		// our triangle start from lower left corner and go counterclockwise

		Vector3d e1 = t->v1 - t->v0;
		Vector3d e2 = t->v2 - t->v0;
		Vector3d nDir = crossProduct(e1, e2); // note the order!
		nDir = normalized(nDir);

		double deltaU1 = t->uv1.x - t->uv0.x;
		double deltaV1 = t->uv1.y - t->uv0.y;

		double deltaU2 = t->uv2.x - t->uv0.x;
		double deltaV2 = t->uv2.y - t->uv0.y;

		double coef = 1 / (-deltaU1 * deltaV2 + deltaV1 * deltaU2);

		Vector3d T = coef * (-deltaV2 * e1 + deltaV1 * e2);
		Vector3d B = coef * (-deltaU2 * e1 + deltaU1 * e2);
		T = normalized(T);
		B = normalized(B);

		Vector3d res;
		res.x = T.x * color.x + B.x * color.y + nDir.x * color.z;
		res.y = T.y * color.x + B.y * color.y + nDir.y * color.z;
		res.z = T.z * color.x + B.z * color.y + nDir.z * color.z;

		inter.nDir = normalized(res);
		break;
	}

	case SPEHRE: {
		Vector3d nDir = inter.nDir;
		Vector3d T = Vector3d(-nDir.y / sqrt(nDir.x * nDir.x + nDir.y * nDir.y),
			nDir.x / sqrt(nDir.x * nDir.x + nDir.y * nDir.y), 0);

		Vector3d B = crossProduct(nDir, T);

		Vector3d res;
		res.x = T.x * color.x + B.x * color.y + nDir.x * color.z;
		res.y = T.y * color.x + B.y * color.y + nDir.y * color.z;
		res.z = T.z * color.x + B.z * color.y + nDir.z * color.z;

		inter.nDir = normalized(res);
		break;
	}

	default:
		break;
	}

}

void textureModify(Intersection& inter, PPMGenerator* g) {
	// DIFFUSE
	if (!FLOAT_EQUAL(-1.0, inter.diffuseIndex)) {
		if (g->diffuseMaps.size() <= inter.diffuseIndex) {
			std::cout <<
				"\ninter.diffuseIndex is greater than diffuseTexuture.size()\nImport texture files in config.txt \n";
			exit(1);
		}
		inter.mtlcolor.diffuse = g->diffuseMaps.at(inter.diffuseIndex)
			->getRGBat(inter.textPos.x, inter.textPos.y);
	}
	// NORMAL
	if (inter.normalMapIndex != -1) {
		changeNormalDir(inter, g);
	}

	// ROUGHNESS
	if (inter.roughnessMapIndex != -1) {
		if (g->roughnessMaps.size() <= inter.diffuseIndex) {
			std::cout <<
				"\ninter.roughnessIndex is greater than roughness_texture.size()\nImport texture files in config.txt \n";
			exit(1);
		}
		inter.mtlcolor.roughness = g->roughnessMaps.at(inter.roughnessMapIndex)
			->getRGBat(inter.textPos.x, inter.textPos.y).x;

	}

	// METALLIC
	if (inter.metallicMapIndex != -1) {
		if (g->metallicMaps.size() <= inter.metallicMapIndex) {
			std::cout <<
				"\ninter.metallicIndex is greater than metallic_texture.size()\nImport texture files in config.txt \n";
			exit(1);
		}
		inter.mtlcolor.metallic = g->metallicMaps.at(inter.metallicMapIndex)
			->getRGBat(inter.textPos.x, inter.textPos.y).x;
	}
}

/// <summary>
/// shadow ray
/// </summary>
/// <param name="p">inter information</param>
/// <param name="lightPos">light position</param>
/// <returns>if the ray is blocked by non-emissive obj, return true</returns>
bool isShadowRayBlocked(Vector3d orig, Vector3d& lightPos, PPMGenerator* g) {
	Vector3d raydir = normalized(lightPos - orig);
	double distance = (lightPos - orig).norm();

	if (!EXPEDITE) {
		Intersection p_light_inter;
		for (auto& i : g->scene.objList) {
			if (i->mtlcolor.hasEmission()) continue; // do not test with light avatar

			if (i->intersect(orig, raydir, p_light_inter) && p_light_inter.t < distance) {
				return true;
			}
		}
		return false;
	}
	else { // BVH intersection test
		return hasIntersection(g->scene.BVHaccelerator->getNode(), orig, raydir, distance);
	}
}

double getLightPdf(Intersection& inter, PPMGenerator* g) {
	if (!inter.intersected) return 0;

	static std::vector<Object*> lightList;
	static double totalArea = 0;
	static bool firstimeCall = true;
	int size = lightList.size();

	// if first time call it, put all the emissive object into lightList

	if (firstimeCall) {
		// 3/2/2024: need lock
		sampleLight_mutex.lock();
		for (auto& i : g->scene.objList) {
			if (!firstimeCall)
				break;

			if (i->mtlcolor.hasEmission()) {
				lightList.emplace_back(i.get());
				// without lock, sometimes problem on i->getArea(), maybe due to unique_ptr
				totalArea += i->getArea();
			}
		}
		firstimeCall = false;
		sampleLight_mutex.unlock();
	}

	size = lightList.size();
	// if there's no light
	if (size == 0) {
		return 0;
	}
	if (!inter.obj->mtlcolor.hasEmission()) return 0;

	// independent event p(a&&b) == p(a) *  p(b)
	double area = inter.obj->getArea();
	return  1 / (size * area);
}

// sample all the emissive object to get one point on their surface,
// update the intersection, and the pdf to sample it
// pdf of that inter is, 1/area of THE object surface area
void sampleLight(Intersection& inter, double& pdf, PPMGenerator* g) {
	static std::vector<Object*> lightList;
	static double totalArea = 0;
	static bool firstimeCall = true;
	int size = lightList.size();

	// if first time call it, put all the emissive object into lightList

	if (firstimeCall) {
		// 3/2/2024: need lock
		sampleLight_mutex.lock();
		for (auto& i : g->scene.objList) {
			if (!firstimeCall)
				break;

			if (i->mtlcolor.hasEmission()) {
				lightList.emplace_back(i.get());
				totalArea += i->getArea();	// without lock, sometimes problem on i->getArea(), maybe due to unique_ptr
			}
		}
		firstimeCall = false;
		sampleLight_mutex.unlock();
	}

	size = lightList.size();
	// if there's no light
	if (size == 0) {
		inter.intersected = false;
		pdf = 0;
		return;
	}

	// 3/3/2024 need a better sample method
	int index = (int)(getRandomFloat() * (size - 1) + 0.4999);
	if (size == 1) index = 0;

	Object* lightObject = lightList.at(index);

	lightObject->samplePoint(inter, pdf);
	// independent event p(a&&b) == p(a) *  p(b)
	pdf = (1.0 / (size * lightObject->getArea()));
}