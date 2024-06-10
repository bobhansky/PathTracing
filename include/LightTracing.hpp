#pragma once
// Integrator: light
#include "IIntegrator.hpp"

#define CHECK 1		// checking world to raster

// light tracing / particle tracing
class LightTracing : public IIntegrator {
public:
	LightTracing(PPMGenerator* g, IIntersectStrategy* inters) {
		this->g = g;
		this->interStrategy = inters;
	}

	bool sampleLightDir(Vector3d& N, double& dirPdf, Vector3d& sampledRes) {
		// cos-weighted
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
			false;

		sampledRes = res;
		return true;
	}

	virtual void integrate(PPMGenerator* g) {
		Camera& cam = g->cam;

#if CHECK
		Vector3d u = crossProduct(cam.fwdDir, cam.upDir);
		u = normalized(u);
		Vector3d v = crossProduct(u, cam.fwdDir);
		v = normalized(v);

		double d = 1;

		// tan(hfov/2) =  nearplane.width/2 : d
		// h/2 = tan(hfov/2) * d 
		// here height_half and width_half are the h and w of the nearplane in the world space
		double width_half = abs(tan(degree2Radians(cam.hfov / 2.0)) * d);
		double aspect_ratio = cam.width / (double)cam.height;
		double height_half = width_half / aspect_ratio;


		// we sample the center of the pixel 
		// so we need to add offset to the center later
		Vector3d n = normalized(g->viewdir);
		Vector3d eyePos = cam.position;
		Vector3d ul = eyePos + d * n - width_half * u + height_half * v;
		Vector3d ur = eyePos + d * n + width_half * u + height_half * v;
		Vector3d ll = eyePos + d * n - width_half * u - height_half * v;
		Vector3d lr = eyePos + d * n + width_half * u - height_half * v;


		Vector3d delta_h = Vector3d(0, 0, 0);	// delta horizontal
		if (g->width != 1) delta_h = (ur - ul) / (g->width - 1);
		Vector3d delta_v = Vector3d(0, 0, 0);	// delta vertical
		if (g->height != 1) delta_v = (ll - ul) / (g->height - 1);
		Vector3d c_off_h = (ur - ul) / (double)(g->width * 2);	// center horizontal offset
		Vector3d c_off_v = (ll - ul) / (double)(g->height * 2); // vertical

		for (int y = 0; y < g->height; y++) {
			Vector3d v_off = y * delta_v;
			//PRINT = false;
			for (int x = 0; x < g->width; x++) {
				if (x == 512 && y == 380) {
					PRINT = true;
				}

				Vector3d h_off = x * delta_h;
				Vector3d pixelPos = ul + h_off + v_off + c_off_h + c_off_v;		// pixel center position in world space
				Vector3d rayDir;
				Vector3d eyeLocation;

				rayDir = normalized(pixelPos - eyePos);
				eyeLocation = eyePos;

				// check if worldPos to pixel works 
				Intersection inter = getIntersection(g->scene.BVHaccelerator->getNode(), eyePos, rayDir);
				if (!inter.intersected)
					continue;
				Vector3d ret = inter.mtlcolor.diffuse;
				Vector3d interPos = inter.pos;
				int i = cam.worldPos2PixelIndex(interPos);
				if (!cam.FrameBuffer.setRGB(i, ret)) {
					//std::cout<<"cam set rgb out of bound\n";
				}
			}
			showProgress((double)y / g->height);
		}
#else
		// sample light
		Intersection vertexInter;
		double pickpdf;
		sampleLight(vertexInter, pickpdf,g);
		// sample direction

#endif

	}
};