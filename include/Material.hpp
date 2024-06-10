#pragma once

#include "global.hpp"
#include "Vector.hpp"
#include <tuple>



enum MaterialType {
	LAMBERTIAN,	// cosine weighted
	PERFECT_REFLECTIVE,
	PERFECT_REFRACTIVE,
	MICROFACET_R,	// MICROFACET_REFLECTIVE, Cook Torrance Microfacet model  with GGX dist
	MICROFACET_T,	// MICROFACET_TRANSMISSIVE
	UNLIT
};


class Material {
public:
	Vector3d diffuse = Vector3d(0.9, 0.9, 0.9);
	Vector3d specular = Vector3d(1.0);
	Vector3d emission = Vector3d(0.0);
	MaterialType mType = LAMBERTIAN;

	double alpha = 1;		// opacity		if alpha = 1. then no refraction 
	double eta = 1;			// index of refraction

	double roughness = 1;	// width parameter  alpha_g
	double metallic = 0;


	// copy assignment operator
	Material& operator= (const Material& other) {
		if (this != &other) {
			this->diffuse.x = other.diffuse.x;
			this->diffuse.y = other.diffuse.y;
			this->diffuse.z = other.diffuse.z;
			this->specular.x = other.specular.x;
			this->specular.y = other.specular.y;
			this->specular.z = other.specular.z;
			this->emission = other.emission;
			this->mType = other.mType;

			this->alpha = other.alpha;
			this->eta = other.eta;

			this->roughness = other.roughness;
			this->metallic = other.metallic;
		}
		return *this;
	}

	bool hasEmission() {
		return emission.x || emission.y || emission.z;
	}

	// BxDF, return vec3f of elements within [0,1] 
	// wi, wo: origin at inter.pos center, pointing outward
	// wi: incident ray
	// wo: view dir
	Vector3d BxDF(const Vector3d& wi, const Vector3d& wo, const Vector3d& N, double eta_scene, bool TIR = false) {
		switch (mType) {
		case LAMBERTIAN: {
			double cos_theta = wo.dot(N);
			// account for reflection contribution only
			if (cos_theta > 0.0) {
				return diffuse / M_PI;
			}
			else return Vector3d(0.0);
			break;
		}

		case MICROFACET_R: {
			// Cook-Torrance Model

			Vector3d h = normalized(wi + wo);
			double costheta = h.dot(wo);

			Vector3d F0(0.04f);	// should be 0.04			
			F0 = lerp(F0, this->diffuse, this->metallic);
			Vector3d F = fresnelSchlick(costheta, F0);	// learnopgl https://learnopengl.com/PBR/Theory
			// double F = fresnel(wo, h, eta_scene, this->eta);
			double D = D_ndf(h, N, roughness);
			//double G = GeometrySmith(N, wi, wo, (roughness + 1) * (roughness + 1) / 8);	// learnopgl
			double G = G_smf(wi, wo, N, roughness, h);
			Vector3d fr = (F * G * D) / (4 * wi.dot(N) * wo.dot(N));	// originaly double fr

			//fr = clamp(0, 1, fr);
			Vector3d diffuse_term = (1.0 - F) * (diffuse / M_PI);
			Vector3d ref_term = fr;
			// return ref_term;	// used for MIS testing
			return diffuse_term + ref_term;
		}
		case MICROFACET_T: {
			double eta_i = eta_scene;
			double eta_t = eta;
			Vector3d interN = N;
			if (wo.dot(N) < 0) {
				interN = -N;
				std::swap(eta_i, eta_t);
			}
			double F = fresnel(wo, interN, eta_i, eta_t);

			// if wi is reflection dir
			if (wi.dot(interN) >= 0) {
				Vector3d h = normalized(wo + wi);
				double cosTheta = h.dot(wo);
				cosTheta = abs(cosTheta);
				double F = fresnel(wo, h, eta_i, eta_t);
				if (TIR) 
					F = 1.0;
				double D = D_ndf(h, interN, roughness);
				double G = G_smf(wi, wo, interN, roughness, h);
				Vector3d fr = (F * G * D) / (4 * wi.dot(interN) * wo.dot(interN));	// originaly double fr
				return fr;
			}
			// wi is refraction dir
			else {	// need h and jacobian 
				Vector3d h = -normalized(eta_i * wo + eta_t * wi);
				if (h.dot(interN) < 0) h = -h;
				double cos_ih = wi.dot(h), cos_oh = wo.dot(h),
					cos_in = wi.dot(interN), cos_on = wo.dot(interN);
				double F = fresnel(wo, h, eta_i, eta_t);
				double D = D_ndf(h, interN, roughness);
				double G = G_smf(wi, wo, interN, roughness, h);
				double numerator = abs(cos_ih) * abs(cos_oh) * eta_t * eta_t * (1 - F) * G * D;
				double denominator = abs(cos_in) * abs(cos_on) * std::pow(eta_i * cos_ih + eta_t * cos_oh, 2);
				return numerator / denominator;
			}
		}

		case PERFECT_REFLECTIVE:{
			if (FLOAT_EQUAL(normalized(wi + wo).dot(N), 1.0))
				// https://www.youtube.com/watch?v=sg2xdcB8M3c
				return 1 /abs(N.dot(wi));
			return 0;
			break;
		}

		case PERFECT_REFRACTIVE: {
			// https://www.youtube.com/watch?v=sg2xdcB8M3c
			// https://www.pbr-book.org/3ed-2018/Reflection_Models/Specular_Reflection_and_Transmission#fragment-BxDFDeclarations-7
			Vector3d refDir = normalized(getReflectionDir(wo, N));
			double eta_i = eta_scene;
			double eta_t = this->eta;
			double F;
			Vector3d interN = N;
			if (wo.dot(N) < 0) {
				interN = -N;
				std::swap(eta_i, eta_t);
			}
			F = fresnel(wo, interN, eta_i, eta_t);
			Vector3d transDir = normalized(getRefractionDir(wo, interN, eta_i, eta_t));

			interN = interN.dot(wi) < 0 ? -interN : interN;

			if (TIR) 
				return 1 / interN.dot(wi);
			if (FLOAT_EQUAL(wi.dot(refDir), 1.0))
				return F * 1 / interN.dot(wi);
			else if(FLOAT_EQUAL(wi.dot(transDir), 1.0)) 
				// 5/13/2024  this term make some faces bright? idk if it's correct
				return /*(eta_t * eta_t) / (eta_i * eta_i) * */(1 - F) * 1 / interN.dot(wi);
			
			return Vector3d(0.0);
			break;			
		}

		default:
			return Vector3d(0.0);
		}
	}



	// sample a direction on the hemisphere, changes the value of "sampledRes"
	// wi: incident dir, pointing outward
	// when passed in, eta_i is always eta_world
	// returns: first bool for sample success, true for succeed
	//			second bool for special event happening, 1 for happened
	std::tuple<bool, bool> sampleDirection(const Vector3d& wi, const Vector3d& N, Vector3d& sampledRes, double eta_i = 0.0) {
		switch (mType)
		{
		case MICROFACET_R: {
			if (wi.dot(N) <= 0.0f) 
				return { false, false };		// crucial
			
			// https://zhuanlan.zhihu.com/p/78146875
			// https://agraphicsguynotes.com/posts/sample_microfacet_brdf/
			double r0 = getRandomFloat();
			double r1 = getRandomFloat();
			double alhpa = roughness * roughness;
			alpha = std::max(alpha, 1e-3);
			double a2 = alhpa * alpha;
			
			double phi = 2 * M_PI * r1;
			double theta = std::acos(sqrt((1 - r0) / (r0 * (a2 - 1) + 1)));

			double r = std::sin(theta);
			Vector3d h = normalized(Vector3d(r * std::cos(phi), r * std::sin(phi), std::cos(theta)));
			Vector3d res = getReflectionDir(wi, SphereLocal2world(N, h));
			res = normalized(res);
			if (res.dot(N) <= 0)	// actually handled in shadow masking term, but only for specular term
				return {false, false};

			sampledRes = res;
			return { true, false };
			break;
		}
		case MICROFACET_T: {
			double r0 = getRandomFloat();
			double r1 = getRandomFloat();
			double a = roughness * roughness;
			a = std::max(a, 1e-3);
			double a2 = a * a;

			double phi = 2 * M_PI * r1;
			double theta = std::acos(sqrt((1 - r0) / (r0 * (a2 - 1) + 1)));

			double r = std::sin(theta);
			Vector3d h = normalized(Vector3d(r * std::cos(phi), r * std::sin(phi), std::cos(theta)));
			double eta_t = eta;
			Vector3d interN = N;
			if (wi.dot(N) < 0) {
				std::swap(eta_i, eta_t);
				interN = -interN;
			}
			h = SphereLocal2world(interN, h);
			
			Vector3d res = getRefractionDir(wi, h, eta_i, eta_t);
			// if TIR
			if (res.norm2() == 0) {
				return { true, true };	// outside will handle sampleDir
			}

			double F = fresnel(wi, h, eta_i, eta_t);
			// choose refl or trans
			if (getRandomFloat() < F) {	
				sampledRes = getReflectionDir(wi, h);
			}
			else 
				sampledRes = res;

			return { true, false };
			break;
		}

		case LAMBERTIAN: {
			if (wi.dot(N) <= 0.0f)
				return {false, false};
			// cosine weighted
			// **** inverse transformation sampling
			// pbrt 13.6.1  *important
			// https://pbr-book.org/3ed-2018/Monte_Carlo_Integration/2D_Sampling_with_Multidimensional_Transformations
			// https://pbr-book.org/3ed-2018/Monte_Carlo_Integration/Transforming_between_Distributions
			// 
			// 3.6 Approximating Distributions
			// https://raytracing.github.io/books/RayTracingTheRestOfYourLife.html#generatingrandomdirections/uniformsamplingahemisphere
			// https://www.youtube.com/watch?v=rnBbYsysPaU&t=1s

			// 1. generate a random direction in sphere coordinate
			// 2. convert it to world corrdinate
			double r1 = getRandomFloat();
			double r2 = getRandomFloat();
			double cosTheta = sqrt(r1);
			double phi = 2 * M_PI * r2;

			Vector3d dir;
			double sinTheta = sqrt(std::max(0.0, 1 - r1));
			dir.x = cos(phi) * sinTheta;
			dir.y = sin(phi) * sinTheta;
			dir.z = cosTheta;

			dir = normalized(dir);

			Vector3d res = SphereLocal2world(N, dir);
			if (normalized(res).dot(N) < 0)
				return { false, false };

			sampledRes = res;
			return { true,false };
			break;
		}

		case PERFECT_REFLECTIVE: {
			sampledRes = getReflectionDir(wi, N);
			return { true,false };
			break;
		}
		case PERFECT_REFRACTIVE: {
			double eta_t = eta;
			Vector3d interN = N;
			if (wi.dot(N) < 0) {
				std::swap(eta_i, eta_t);
				interN = -interN;
			}

			Vector3d res = getRefractionDir(wi, interN, eta_i, eta_t);
			// if TIR
			if (res.norm2() == 0) {
				return { true, true };	// outside will handle sampleDir
			}

			double F = fresnel(wi, interN, eta_i, eta_t);
			if (getRandomFloat() < F) {
				sampledRes = getReflectionDir(wi, interN);
			}
			else sampledRes = res;

			return { true,false };
			break;
		}
		default: {
			return { false,false };
			break;
		}
		}
	}


	// wo: -camera dir   wi: sampled dir
	// when passed in, eta_i is always eta_world, eta_t is always ior of inter.material
	double pdf(const Vector3d& wi, const Vector3d& wo, const Vector3d& N, double eta_i = 0.0, double eta_t = 0.0) const {
		switch (mType)
		{
		case LAMBERTIAN: {
			// uniform sample probability 1 / (2 * PI)
			if (wi.dot(N) > 0.0)
				return wi.dot(N) / M_PI;  // cosine weighted pdf https://ameye.dev/notes/sampling-the-hemisphere/
			else
				return 0.0;

			break;
		}
		case MICROFACET_R: {
			// corresponds to normal distribution function D
			// https://www.tobias-franke.eu/log/2014/03/30/notes_on_importance_sampling.html
			Vector3d h = normalized(wo + wi);
			double cosTheta = N.dot(h);
			cosTheta = std::max(cosTheta, 0.0);

			return D_ndf(h, N, roughness) * cosTheta / (4.0 * wo.dot(h));
			break;
		}
		case MICROFACET_T: {
			Vector3d interN = N;
			if (wo.dot(N) < 0) {
				interN = -N;
				std::swap(eta_i, eta_t);
			}
			double F = fresnel(wo, interN, eta_i, eta_t);

			// if wi is reflection dir
			if (wi.dot(interN) >= 0) {
				Vector3d h = normalized(wo + wi);
				double cosTheta = interN.dot(h);
				cosTheta = abs(cosTheta);
				return F *  D_ndf(h, interN, roughness) * cosTheta / (4.0 * wo.dot(h));
			}
			// wi is refraction dir
			else {	// need h and jacobian 
				Vector3d h = -normalized(eta_i * wo + eta_t * wi);
				double cosTheta = interN.dot(h);
				if (cosTheta < 0) {
					h = -h;
					cosTheta = abs(cosTheta);
				}
				double denominatorSqrt = eta_i * wi.dot(h) + eta_t * wo.dot(h);
				double jacobian = (eta_t * eta_t * abs(wo.dot(h)))/ (denominatorSqrt * denominatorSqrt);
				return (1 - F) *  D_ndf(h, interN, roughness) * cosTheta * jacobian;
			}
			break;
		}

		case PERFECT_REFLECTIVE: {
			if(FLOAT_EQUAL(normalized(wi+wo).dot(N), 1.0))
				return 1;
			 return 0;
			break;
		}

		case PERFECT_REFRACTIVE: {
			Vector3d refDir = normalized(getReflectionDir(wo, N));
			Vector3d nDir = N;
			if (wo.dot(nDir) < 0) {
				std::swap(eta_i, eta_t);
				nDir = -N;
			}
			Vector3d transDir = normalized(getRefractionDir(wo, nDir, eta_i, eta_t));

			double F = fresnel(wo, nDir, eta_i, eta_t);

			// check which direction is sampled 
			if (FLOAT_EQUAL(wi.dot(refDir), 1.0))
				return F;
			else if (FLOAT_EQUAL(wi.dot(transDir), 1.0))
				return 1 - F;

			return 0;
			break;
		}

		default:
			return 1;
			break;
		}
	}
};