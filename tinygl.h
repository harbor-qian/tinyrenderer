#pragma once

#include "model.h"
#include "geometry.h"
#include <windows.h>

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green = TGAColor(0, 255, 0, 255);

extern const int width;
extern const int height;

extern Matrix ModelView;
extern Matrix Projection;
extern Matrix Viewport;

struct IShader {
    virtual ~IShader();
    virtual Vec4f vertex(int iface, int nthvert) = 0;
    virtual bool fragment(Vec3f bar, TGAColor& color) = 0;
};

void viewport(int x, int y, int w, int h);
void projection(float fov, float aspect, float n, float f);
void lookat(Vec3f eye, Vec3f center, Vec3f up);

void line(int x0, int y0, int x1, int y1, TGAImage& image, TGAColor color);

void triangle_s(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage& image, TGAColor color);

Vec3f barycentric(Vec3f* pts, Vec3f P);

float EdgeFunc(const Vec3f& p0, const Vec3f& p1, const Vec3f& p2);

//void triangle(Vec4f* pts, IShader& shader, TGAImage& image, TGAImage& zbuffer);
void triangle(Vec4f* pts, IShader& shader, TGAImage& image, float* zbuffer);


