#pragma once
#include "global.hpp"
#include "Vector.hpp"

#include <vector>


class Texture {
public:
	std::string name;
	int width = 0;
	int height = 0;
	std::vector<Vector3d> rgb;

	// by u v
	Vector3d getRGBat(double u, double v) {
		if (width == 0 && height == 0) {
			return Vector3d();
		}

		int x = u * width;
		int y = v * height;
		// double precision might lead to out of bounds 
		// so clamp it
		int index = y * width + x;
		if (index < 0) index = 0;
		if (index >= rgb.size()) index = rgb.size()-1;
		return getEleIn(rgb, index);
	}

	bool setRGB(int x, int y, Vector3d& RGB) {
		if (x < 0 || x >= width || y < 0 || y >= height)
			return false;

		rgb[y + x * width] = RGB;
		return true;
	}

	bool setRGB(int index, Vector3d& RGB) {
		if (index < 0 || index >= width * height)
			return false;

		rgb[index] = RGB;
		return true;
	}
};