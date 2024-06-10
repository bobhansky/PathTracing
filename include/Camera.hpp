#include "Vector.hpp"
#include "Texture.hpp"

// matrix is copied from SmallVCM

class Camera {
public:
	Camera() {
	}

	// initialize all the data
	void initialize(Vector3d bkgcolor) {
		fwdDir = normalized(fwdDir);
		rightDir = crossProduct(fwdDir, upDir);
		rightDir = normalized(rightDir);
		upDir = crossProduct(rightDir, fwdDir);
		upDir = normalized(upDir);


		// projection position onto view basis
		Vector3d pos(
			rightDir.dot(position),
			upDir.dot(position),
			(-fwdDir).dot(position)
		);

		FrameBuffer.rgb.resize(width * height);
		FrameBuffer.rgb.assign(width * height, Vector3d(bkgcolor.x, bkgcolor.y, bkgcolor.z));	//  set to bkgcolor
		FrameBuffer.width = width;
		FrameBuffer.height = height;

		world2Cam.setRow(0, rightDir, -pos.x);
		world2Cam.setRow(1, upDir, -pos.y);
		world2Cam.setRow(2, -fwdDir, -pos.z);
		world2Cam.setRow(3, Vector3d(0.0), 1.0);

		perspective = getPerspectiveMatrix(hfov, 0.1, 10000.0, (double)width/height);
		world2ndc = perspective * world2Cam;
		world2Raster = Mat4d::getTranslate(Vector3d(1.0, 1.0, 0)) * world2ndc;
		world2Raster = Mat4d::getScale(Vector3d(width * 0.5, height * 0.5, 0)) * world2Raster;
	}

	// return -1 if can't
	int raster2pxlIndex(Vector4d& raster) {
		int x = (int)raster.x;
		int y = (int)raster.y;
		if (x < 0 || x >= width || y < 0 || y >= height)
			return -1;

		return x + width * y;
	}

	int worldPos2PixelIndex(const Vector3d& pos) {
		// Vector3f raster1 = world2Raster.transformPoint(pos); 
		Vector4d raster  = world2Raster * Vector4d(pos.x , pos.y, pos.z, 1);
		raster.normalizeW();
		// ray generation is offset by 0.5 pixel, now offset it back
		raster.x -= 0.5;
		raster.y -= 0.5;
		
		// debug
		//Vector4f pCam = world2Cam * Vector4f(pos.x, pos.y, pos.z, 1);
		//pCam.normalizeW();
		//Vector4f pNDC = perspective * pCam;
		//pNDC.normalizeW();
		//Vector4f pR = Mat4f::getTranslate(Vector3f(1.0, 1.0, 0)) * pNDC;
		//pR.normalizeW();
		//pR = Mat4f::getScale(Vector3f(width * 0.5, height * 0.5, 0)) * pR;

		return raster2pxlIndex(raster);
	}

public:
	Vector3d fwdDir;
	Vector3d upDir;
	Vector3d rightDir;
	Vector3d position;
	Texture FrameBuffer;
	int width;
	int height;
	int hfov;

	Mat4d world2Cam;
	Mat4d perspective;
	Mat4d world2ndc;
	Mat4d world2Raster;
};