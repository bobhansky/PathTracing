#pragma once

#include <math.h>
#include <iostream>




class Vector3i {
public:
	int x;
	int y;
	int z;

	Vector3i(int xval, int yval, int zval) : x(xval), y(yval),z(zval) {}

	Vector3i() {
		x = 0;
		y = 0;
		z = 0;
	}
};

class Vector4d {
public:
	double x;
	double y;
	double z;
	double w;

	Vector4d(double xval, double yval, double zval, double wval) : x(xval), y(yval), z(zval), w(wval) {}

	Vector4d() {
		x = 0.0;
		y = 0.0;
		z = 0.0;
		w = 0.0;
	}

	void normalizeW() {

		x = x / w;
		y = y / w;
		z = z / w;
		w = w / w;
	}
};

class Vector2d {
public:
	double x;
	double y;
	Vector2d(double xval, double yval) : x(xval), y(yval) {}
	Vector2d() {
		x = -1.0;
		y = -1.0;
	}

	Vector2d operator*(const double& c) {
		return Vector2d(x * c, y * c);
	}
	Vector2d operator+(const Vector2d& v) const
	{
		return Vector2d(x + v.x, y + v.y);
	}
};



class Vector3d {
public:
	double x;
	double y;
	double z;

	Vector3d(double xval, double yval, double zval) : x(xval), y(yval), z(zval) {}

	Vector3d() {
		x = 0.0;
		y = 0.0;
		z = 0.0;
	}

	Vector3d(double i) {		// implicit 	Vector3f a = 1;
		x = i;
		y = i;
		z = i;
	}
	double get(int i) const {
		switch (i)
		{
		case 0: {
			return x;
			break;
		}
		case 1: {
			return y;
			break;
		}
		case 2: {
			return z;
			break;
		}
		default:
			return -1;
			break;
		}
	}

	void print() {
		std::cout << "x: " << x << ", y: " << y << ", z: " << z << std::endl;
	}

	std::string toString() const {
		return "(" + std::to_string(x)+ ", "+ std::to_string(y) +", "+ std::to_string(z) + ")";
	}

	// ************************* vector operations *************************

	operator std::string() const {
		return "(" + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z);
	}
	
	// const double&		const lvalue reference, can point to both lvalue and rvalue
	Vector3d operator*(const double& c) {
		return Vector3d(x * c, y * c, z * c);
	}

	Vector3d operator/(const double& c) const{
		return Vector3d(x / c, y / c, z / c);
	}

	// element-wise product
	Vector3d operator*(const Vector3d& v) const
	{
		return Vector3d(x * v.x, y * v.y, z * v.z);
	}

	Vector3d operator/(const Vector3d& v) const {
		return Vector3d(x / v.x, y / v.y, z / v.z);
	}

	Vector3d operator-(const Vector3d& v) const
	{
		return Vector3d(x - v.x, y - v.y, z - v.z);
	}
	Vector3d operator+(const Vector3d& v) const
	{
		return Vector3d(x + v.x, y + v.y, z + v.z);
	}
	void operator+=(const Vector3d& v)
	{
		x = x + v.x;
		y = y + v.y;
		z = z + v.z;
	}
	void operator-=(const Vector3d& v)
	{
		x = x - v.x;
		y = y - v.y;
		z = z - v.z;
	}
	Vector3d operator -() const 
	{ 
		return Vector3d(-x, -y, -z); 
	}

	// copy assignment operator
	Vector3d& operator= (const Vector3d& other) {
		this->x = other.x;
		this->y = other.y;
		this->z = other.z;

		return *this;
	}
	
	Vector3d& operator= (const double& scaler) {
		this->x = scaler;
		this->y = scaler;
		this->z = scaler;

		return *this;
	}

	// dot product
	double dot(const Vector3d& v) const {
		return x * v.x + y * v.y + z * v.z;
	}

	friend Vector3d operator*(double c, const Vector3d& v) {
		return Vector3d(v.x * c, v.y * c, v.z * c);
	}

	friend Vector3d operator-(double c, const Vector3d& v) {
		return Vector3d(c - v.x, c - v.y, c - v.z);
	}


	// ************************* vector operations ends ***********************

	// get norm of this vector
	double norm() {
		return sqrt(x * x + y * y + z * z);
	}

	double norm2() {
		return x * x + y * y + z * z;
	}

};

// return the normalized version of vector v
inline Vector3d normalized(const Vector3d& v) {
	double mag = sqrt((v.x * v.x + v.y * v.y + v.z * v.z));
	if (mag > 0) {
		double mag_inv = 1 / mag;		// for efficiency we times instead of divide 3 times
		return Vector3d(v.x * mag_inv, v.y * mag_inv, v.z * mag_inv);
	}
	return v;
}

// return the crossProduct of v1 and v2
Vector3d crossProduct(const Vector3d& v1, const Vector3d& v2) {
	return Vector3d(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}


class Mat4d {
public:
	Mat4d() {
		for (int i = 0; i < 16; i++)
			ele[i] = 0;
	}

	Mat4d(double f1, double f2, double f3, double f4, double f5,
		double f6, double f7, double f8, double f9, double f10,
		double f11, double f12, double f13, double f14, double f15, double f16 ) {
		ele[0] = f1; ele[1] = f2; ele[2] = f3; ele[3] = f4;
		ele[4] = f5; ele[5] = f6; ele[6] = f7; ele[7] = f8;
		ele[8] = f9; ele[9] = f10; ele[10] = f11; ele[11] = f12;
		ele[12] = f13; ele[13] = f14; ele[14] = f15; ele[15] = f16;
	}

	double get(int r, int c) const {
		return ele[c + r * 4];
	}

	void set(int r, int c, double val) {
		ele[c + r * 4] = val;
	}
	void setRow(int r, const Vector3d& vec, double val) {
		ele[0 + r * 4] = vec.x;
		ele[1 + r * 4] = vec.y;
		ele[2 + r * 4] = vec.z;
		ele[3 + r * 4] = val;
	}

	void print() {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
			std::cout << get(i, j) << " ";
			}
			std::cout << "\n";
		}
	}


	

	static Mat4d getTranslate(const Vector3d& v){
		Mat4d r;
		r.set(3, 3, 1);
		r.set(0, 3, v.x);
		r.set(1, 3, v.y);
		r.set(2, 3, v.z);

		r.set(0, 0, 1); r.set(1, 1, 1); r.set(2, 2, 1); r.set(3, 3, 1);
		return r;
	}

	static Mat4d getScale(const Vector3d& v){
		Mat4d r;
		r.set(3, 3, 1);
		r.set(0, 0, v.x);
		r.set(1, 1, v.y);
		r.set(2, 2, v.z);
		return r;
	}

	Vector4d operator*(const Vector4d& v) {
		Vector4d res;
		res.x = v.x * ele[0] + v.y * ele[1] + v.z * ele[2] + v.w * ele[3];
		res.y = v.x * ele[4] + v.y * ele[5] + v.z * ele[6] + v.w * ele[7];
		res.z = v.x * ele[8] + v.y * ele[9] + v.z * ele[10] + v.w * ele[11];
		res.w = v.x * ele[12] + v.y * ele[13] + v.z * ele[14] + v.w * ele[15];
		return res;
	}

	// copy from smallVCM
	Vector3d transformPoint(const Vector3d& p) {
		// get calculated W,vector.w, the forth dimension
		double w = get(3, 3);

		for (int c = 0; c < 3; c++)
			w += get(3, c) * p.get(c);

		// normalization factor
		const double invW = 1.0 / w;

		Vector3d res(0);
		// translate value
		res.x = get(0, 3);
		// dot product
		for (int c = 0; c < 3; c++)
			res.x += p.get(c) * get(0, c);
		res.x *= invW;

		res.y = get(1, 3);
		for (int c = 0; c < 3; c++)
			res.y += p.get(c) * get(1, c);
		res.y *= invW;

		res.z = get(2, 3);
		for (int c = 0; c < 3; c++)
			res.z += p.get(c) * get(2, c);
		res.z *= invW;

		return res;
	}


public:
	double ele[16];
};

// copy from SmallVCM
Mat4d operator*(const Mat4d& left, const Mat4d& right)
{
	Mat4d res;
	for (int row = 0; row < 4; row++) {
		for (int col = 0; col < 4; col++) {
			double r = 0;
			for (int i = 0; i < 4; i++)
				r += left.get(row, i) * right.get(i, col);
			res.set(row, col, r);
		}
	}
	return res;
}

// copy from SmallVCM
Mat4d getPerspectiveMatrix(double aFov, double aNear, double aFar, double aspect_ratio) {
	// Camera points towards -z.  0 < near < far.
	// Matrix maps z range [-near, -far] to [-1, 1], after homogeneous division.
	double f = 1.0 / (std::tan(aFov * 3.1415926535897 / 360.0));
	double d = 1.0 / (aNear - aFar);
	Mat4d r;
	r.ele[0] = f;     r.ele[1] = 0.0;  r.ele[2] = 0.0;					r.ele[3] = 0.0;
	r.ele[4] = 0.0;  r.ele[5] = -f;    r.ele[6] = 0.0;					r.ele[7] = 0.0;
	r.ele[8] = 0.0;  r.ele[9] = 0.0;  r.ele[10] = (aNear + aFar) * d;		r.ele[11] = 2.0 * aNear * aFar * d;
	r.ele[12] = 0.0; r.ele[13] = 0.0; r.ele[14] = -1.0;					r.ele[15] = 0.0;
	return r;


	//Mat4f p2o;
	//p2o.ele[0] = aNear;     p2o.ele[1] = 0.0f;  p2o.ele[2] = 0.0f;					p2o.ele[3] = 0.0f;
	//p2o.ele[4] = 0.0f;  p2o.ele[5] = aNear;     p2o.ele[6] = 0.0f;					p2o.ele[7] = 0.0f;
	//p2o.ele[8] = 0.0f;  p2o.ele[9] = 0.0f;		p2o.ele[10] = (aNear + aFar);			p2o.ele[11] = aNear * aFar;
	//p2o.ele[12] = 0.0f; p2o.ele[13] = 0.0f;		p2o.ele[14] = -1.0f;					p2o.ele[15] = 0.0f;
	//double t = std::tan((aFov / 2) * 3.1415926535897f / 180) * aNear;
	//double b = -t;
	//double r = aspect_ratio * t;
	//double l = -r;
	//Mat4f orth_trans(1, 0, 0, -(r + l) / 2,
	//	0, 1, 0, -(t + b) / 2,
	//	0, 0, 1, -(aNear + aFar) / 2,
	//	0, 0, 0, 1);
	//Mat4f orth_scale(2 / (r - l), 0, 0, 0,
	//	0, 2 / (t - b), 0, 0,
	//	0, 0, 2 / (aNear - aFar), 0,
	//	0, 0, 0, 1);
	//Mat4f orth = orth_scale * orth_trans;
	//Mat4f proj = orth * p2o;
	//return proj;
}