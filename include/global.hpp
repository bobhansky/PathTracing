#pragma once
#include <string>
#include <float.h>
#include <stdexcept>
#include <math.h>
#include <algorithm>
#include <random>
#include <iostream>
#include <stdio.h>

#include "Vector.hpp"


#define M_PI 3.1415926535897
#define EPSILON 0.05		// be picky about it, change it to accommodate object size

// lerp(x,v0,v1) = v0 + x(v1-v0);
// x is the portion
// v0 v1 are the values
double lerp(double v0, double v1, double x) {
	return v0 + x * (v1 - v0);
}

Vector3d lerp(Vector3d& v0, Vector3d& v1, double x) {
	Vector3d res;
	res.x = v0.x + x * (v1.x - v0.x);
	res.y = v0.y + x * (v1.y - v0.y);
	res.z = v0.z + x * (v1.z - v0.z);

	return res;
}

inline double clamp(const double& lo, const double& hi, const double& v)
{
	return std::max(lo, std::min(hi, v));
}

inline Vector3d clamp(const Vector3d& lo, const Vector3d& hi, Vector3d& v)
{
	Vector3d res;
	res.x = clamp(lo.x, hi.x, v.x);
	res.y = clamp(lo.y, hi.y, v.y);
	res.z = clamp(lo.z, hi.z, v.z);
	return res;
}

inline double rescale(double input, double originMax, double originMin, double targetMax, double targetMin) {
	return targetMin + ((targetMax - targetMin) * (input - originMin)/(originMax - originMin));
}


// check if str is convertable to a positive integer
// if there's a char other than 0 to 9
// then throw exeption
void checkPosInt(std::string& str) {
	for (auto i : str) {
		if (i < 48 || i > 57) {
			throw std::runtime_error(str + ": expect a positive number");
		}
	}
}

// check if str is convertable to a double
void checkFloat(std::string& str) {
	if (str.size() == 0) {
		throw std::runtime_error(str + ": not a valid double number");
	}

	bool dotAppeared = false;
	// for minus 
	if (str.at(0) == '-') {
		if(str.size()==1) throw std::runtime_error(str + ": not a valid double number");

		for (int i = 1; i < str.size(); i++) {
			// -abcdef   a can't be a char other than a number
			if (i == 1 && str[1] < 48 || str[i] > 57) throw std::runtime_error(str + ": not a valid double number");
			// for other index 
			else {
				// dot can only present once && can't be at the end of the str
				if (str[i] == '.' && !dotAppeared && i!=str.size()-1) {
					dotAppeared = true;		// 
				}
				else {	// expect numbers
					if(str[i] < 48 || str[i] > 57) throw std::runtime_error(str + ": not a valid double number");
				}
			}
		}
	}
	// for positive
	else {
		for (int i = 0; i < str.size(); i++) {
			// abcdef   a can't be a char other than a number
			if (i == 0 && str[0] < 48 || str[i] > 57) throw std::runtime_error(str + ": not a valid double number");
			// for other index 
			else {
				// dot can only present once && can't be at the end of the str
				if (str[i] == '.' && !dotAppeared && i != str.size() - 1) {
					dotAppeared = true;		// 
				}
				else {	// expect numbers
					if (str[i] < 48 || str[i] > 57) throw std::runtime_error(str + ": not a valid double number");
				}
			}
		}
	}

}

// convert degree to radians
double degree2Radians(const double& d) {
	return d * M_PI / 180;
}

// check if two double numbers are equal
inline bool FLOAT_EQUAL(const double& x, const double& y) {
	return (abs(x - y) < 0.00001);
}

/// <summary		// ctrl + / 
/// solve
/// A*t^2 + B*t + C = 0
/// </summary>
/// <param name="t1"> first solution  </param>
/// <param name="t2"> second solution	</param>
/// 
/// if there's only one solution, then either t1 or t2 == FLT_MAX
/// if there's no real solution, t1 == t2 == FLT_MAX
void solveQuadratic(double& t1, double& t2, double& A, double& B, double& C) {
	double discriminant = B * B - 4 * A * C;
	// no real solution
	if (discriminant < 0) {
		t1 = FLT_MAX;
		t2 = FLT_MAX;
	}
	// one real solution
	else if (discriminant == 0) {
		t1 = (-B + sqrt(discriminant)) / (2 * A);
		t2 = t1;
	}
	// two real solution
	else {
		t1 = (-B + sqrt(discriminant)) / (2 * A);
		t2 = (-B - sqrt(discriminant)) / (2 * A);
	}
	
	// always make t1 the smaller result
	if (t1 > t2) std::swap(t1, t2);
}



// check if a ele existing within a string vector
bool existIn(std::string& ele, std::vector<std::string>& v) {
	for (auto i : v) {
		if (ele.compare(i) == 0) return true;
	}
	return false;
}



// get a uniformly distributed number in range [0,1)
double getRandomFloat() {
	// see 
	// https://stackoverflow.com/questions/38367976/do-stdrandom-device-and-stdmt19937-follow-an-uniform-distribution

	// an uniformly - distributed random number generator, use it to seed a pseudo-random generator
	static std::random_device dev;
	// a fast pseudo-random number generator, use this to seed a particular distribution
	static std::mt19937 rng(dev());		
	static std::uniform_real_distribution<double> dist(0,1); // distribution in range [0.0, 1.0)

	return dist(rng);	
}

// cout to terminal the progress
void showProgress(double prog) {
	int barWidth = 60;

	int pos = barWidth * prog;
	for (int i = 0; i < barWidth; ++i) {
		if (i < pos) std::cout << "=";
		else if (i == pos) std::cout << ">";
		else std::cout << " ";
	}

	std::cout << int(prog * 100.0) << " %\r";
}


// safely get vec3f element in an vec3f array
Vector3d getEleIn(std::vector<Vector3d>& arr, int index) {
	if (index >= arr.size() || index < 0) {
		throw std::runtime_error(index + " is out of bound: array has size: " + arr.size());
	}
	return arr.at(index);
}

// for uv, texture mapping
Vector2d getEleIn(std::vector<Vector2d>& arr, int index) {
	if (index >= arr.size() || index < 0) {
		throw std::runtime_error(index + " is out of bound: array has size: " + arr.size());
	}
	return arr.at(index);
}

// Schlick approximation used for microfacet model
// https://learnopengl.com/PBR/Theory:
// For conductor surfaces (metals), calculating the base reflectivity with indices of refraction doesn't properly hold 
// and we need to use a different Fresnel equation for conductors altogether. 
Vector3d fresnelSchlick(double cosTheta, const Vector3d& F0)
{
	return F0 + (1.0 - F0) * pow(1.0 - cosTheta, 5.0);
}

// fresnel, get the specular reflection fraction 
double fresnel(const Vector3d& Incident, const Vector3d& normal, const double eta_i, const double eta_t) {
	Vector3d I = normalized(Incident);
	Vector3d N = normalized(normal);
	// there are two possible cases:
	// 1. ray is bouncing at outer surface
	// 2. ray os bouncing at innner surface
	// if ray is bouncing at inner surface (I dot N < 0)
	// reverse N direction
	double cosI_N = I.dot(N);
	if (cosI_N < 0) N = -N;
	
	// The Schlick approximation defines the Fresnel
	// reflectance coefficient using the function :
	// Fr = F0 + (1–F0 )(1–cos(theta_i))^5
	// Schlick approximation: a faster approach to define F0
	double F0 = pow(((eta_t - eta_i) / (eta_t + eta_i)), 2.0);
	double Fr = F0 + (1 - F0) * (pow(1 - (I.dot(N)), 5));

	return Fr;
}

// get the reflection direction (un-normalized) by given incident ray and inter's normal direction 
Vector3d getReflectionDir(const Vector3d& incident, const Vector3d& normal) {
	Vector3d I = normalized(incident);
	Vector3d N = normalized(normal);
	
	return 2 * (N.dot(I)) * N - I;
}

// get the transmittance ray direction  (un-normalized)
Vector3d getRefractionDir(const Vector3d& incident, const Vector3d& normal, double eta_i, double eta_t) {
	// see 03-15.raytracing.pdf page 68
	// additional notes:
	// there are two possible cases:
	// 1. ray is traveling from outside to inside of the obj
	// 2. ray is traveling from inside of the obj to outside
	// since we define object's normal pointing toward the outside,
	// we can check the sign of I dot N to tell which case it is
	Vector3d I = normalized(incident);
	Vector3d N = normalized(normal);
	double cos_theta_i = N.dot(I);
	cos_theta_i = clamp(-1, 1, cos_theta_i);

	if (cos_theta_i < 0) {
		N = -N;
		cos_theta_i = -cos_theta_i;
	}

	double sin_theta_i = sqrt(1 - pow(cos_theta_i, 2));
	double sin_theta_t = (eta_i / eta_t) * sin_theta_i;

	// check total internal reflection case:
	// if it is the case, return vec3f(0), meanning no refraction dir
	if (sin_theta_i > (eta_t / eta_i)) 
		return Vector3d(0);

	double cos_theta_t = sqrt(1 - pow(sin_theta_t, 2));

	return cos_theta_t * (-N) + eta_i / eta_t * (cos_theta_i * N - I);
}

/// <summary>
/// Normal Distribution Function
/// isotropic GGX
/// </summary>
/// <param name="h">: half vector, also the micro normal m </param>
/// <param name="n">: macrosurface normal  </param>
/// <param name="roughness">: width parameter alpha_g </param>
/// <returns>return the area of microfacet with micro normal h in dA</returns>
double D_ndf(const Vector3d& h, const Vector3d& n, double roughness) {
	double alpha = roughness * roughness;
	alpha = std::max(alpha, 1e-3);
	if (n.dot(h) < 0) 
		return 0;
	double cos_nh_2 = (n.dot(h)) * (n.dot(h));
	double sin_nh_2 = 1 - cos_nh_2;
	double sum = alpha * alpha * cos_nh_2 + sin_nh_2;
	if (sum == 0)
		return 1;
	double res = (alpha * alpha) / (M_PI * (sum * sum));

	return res;
}

/// <summary>
/// shadow masking function
/// </summary>
/// <param name="wi">: incident solid angle</param>
/// <param name="wo">: observing/out solid angle</param>
/// <param name="n">: normal</param>
/// <param name="roughness">: width parameter alpha_g </param>
/// <returns> the fraction of unblocked part, [0,1]</returns>
double G_smf(const Vector3d& wi, const Vector3d& wo, const Vector3d& n, double roughness, const Vector3d& h) {
	double alpha = roughness * roughness;
	alpha = std::max(alpha, 1e-3);
	double angle_wi_n = acos(wi.dot(n));
	double angle_wo_n = acos(wo.dot(n));
	// in paper   Microfacet Models for Refraction through Rough Surface
	double G1_wi = ((wi.dot(h)/wi.dot(n)) < 0 ? 0 : 1) * 2 / (1 + sqrt(1 + alpha * alpha * pow(tan(angle_wi_n), 2)));
	double G1_wo = ((wo.dot(h)/wo.dot(n)) < 0 ? 0 : 1) * 2 / (1 + sqrt(1 + alpha * alpha * pow(tan(angle_wo_n), 2)));

	return G1_wi * G1_wo;

	
	// a better G  according to https://zhuanlan.zhihu.com/p/434964126
	//double Ai = (-1 + sqrt(1 + roughness * roughness * pow(tanf(angle_wi_n), 2))) * 0.5;
	//double Ao = (-1 + sqrt(1 + roughness * roughness * pow(tanf(angle_wo_n), 2))) * 0.5;
	//return 1 / (1 + Ai + Ao);
}

// from https://learnopengl.com/PBR/Theory
double GeometrySchlickGGX(double NdotV, double k)
{
	double nom   = NdotV;
	double denom = NdotV * (1.0 - k) + k;
	
    return nom / denom;
}

// V 
double GeometrySmith(Vector3d& N, Vector3d& wi, Vector3d& wo, double k)
{
	double NdotV = std::max(N.dot(wo), 0.0);
	double NdotL = std::max(N.dot(wi), 0.0);
	double ggx1 = GeometrySchlickGGX(NdotV, k);
	double ggx2 = GeometrySchlickGGX(NdotL, k);

	return ggx1 * ggx2;
}

double getMisWeight(double pdf, double otherPdf) {
	// balance heuristic
	//return pdf / (pdf + otherPdf);

	// power heuristic
	return (pdf * pdf) / ((pdf + otherPdf) * (pdf + otherPdf));
}

// avoid self colission when testing ray intersection
void offsetRayOrig(Vector3d& orig, Vector3d interNormal, bool rayIsInside = false) {
	rayIsInside ? orig -= interNormal * EPSILON : orig += interNormal * EPSILON;
}

Vector3d SphereLocal2world(const Vector3d& n, const Vector3d& dir) {
	// https://raytracing.github.io/books/RayTracingTheRestOfYourLife.html#generatingrandomdirections/uniformsamplingahemisphere
	// 8. orthonormal basis
	// change of basis 

	// x y z local coordinates to s t n coordinates
	Vector3d a;
	Vector3d N = normalized(n);	//z
	// construct an Orthonalmal basis
	// randomly choose an a that is not parallel to N
	if (abs(N.x) > 0.9)
		a = { 0.0, 1.0, 0.0 };
	else a = { 1.0, 0.0, 0.0 };
	//Vector3f T = crossProduct(a, N); // y   X cross Y == Z      then S cross T should == N
	//Vector3f S = crossProduct(T, N); // x
	// ******** 
	// 2/21/2024 IMPORTANT
	// 2 unit vectors cross product doens't guarantee to produce unit vec, unless they are orthogonal
	// Vector3f S = crossProduct(N, a);		// reason for wrong result
	Vector3d S = normalized(crossProduct(N, a));
	Vector3d T = crossProduct(N, S);

	return normalized(dir.x * S + dir.y * T + dir.z * N);
}
