﻿#include "tinygl.h"

IShader::~IShader() {}


void viewport(int x, int y, int w, int h) {
    Viewport = Matrix::identity();
    Viewport[0][3] = x + w / 2.f;
    Viewport[1][3] = y + h / 2.f;
    Viewport[2][3] = 255.f / 2.f;
    Viewport[0][0] = w / 2.f;
    Viewport[1][1] = h / 2.f;
    Viewport[2][2] = 255.f / 2.f;
}

void projection(float coeff) {
    Projection = Matrix::identity();
    Projection[3][2] = coeff;
}

void lookat(Vec3f eye, Vec3f center, Vec3f up) {
    Vec3f z = (eye - center).normalize();
    Vec3f x = cross(up, z).normalize();
    Vec3f y = cross(z, x).normalize();
    ModelView = Matrix::identity();
    for (int i = 0; i < 3; i++) {
        ModelView[0][i] = x[i];
        ModelView[1][i] = y[i];
        ModelView[2][i] = z[i];
        ModelView[i][3] = -center[i];
    }
}


void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) {
    bool steep = false;
    if (std::abs(x0-x1)<std::abs(y0-y1)) {
        std::swap(x0, y0);
        std::swap(x1, y1);
        steep = true;
    }
    if (x0>x1) {
        std::swap(x0, x1);
        std::swap(y0, y1);
    }
    int dx = x1 - x0;
    int dy = y1 - y0;

    int y = y0;
    float error = 0;
    int yincr = y1 > y0 ? 1 : -1;
    int derror = std::abs(dy) * 2;

    if (steep) {
        for (int x = x0; x <= x1; x++) {
            image.set(y, x, color);
            error += derror;
            if (error > dx) {
                error -= 2 * dx;
                y += yincr;
            }
        }
    }
    else {
        for (int x = x0; x <= x1; x++){
            image.set(x, y, color);
            error += derror;
            if (error > dx) {
                error -= 2 * dx;
                y += yincr;
            }
        }
    }
}

void triangle_s(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage& image, TGAColor color) {
    // sort the vertices, t0, t1, t2 lower−to−upper (bubblesort yay!) 
    if (t0.y > t1.y) std::swap(t0, t1);
    if (t0.y > t2.y) std::swap(t0, t2);
    if (t1.y > t2.y) std::swap(t1, t2);
    int total_height = t2.y - t0.y;
    for (int y = t0.y; y <= t1.y; y++) {
        int segment_height = t1.y - t0.y + 1;
        float alpha = (float)(y - t0.y) / total_height;
        float beta = (float)(y - t0.y) / segment_height; // be careful with divisions by zero 
        Vec2i A = t0 + (t2 - t0) * alpha;
        Vec2i B = t0 + (t1 - t0) * beta;
        line(A.x, y, B.x, y, image, green);
        //if (A.x > B.x) std::swap(A, B);
        //for (int j = A.x; j <= B.x; j++) {
        //    image.set(j, y, color); // attention, due to int casts t0.y+i != A.y 
        //}
    }
    for (int y = t1.y; y <= t2.y; y++) {
        int segment_height = t2.y - t1.y + 1;
        float alpha = (float)(y - t0.y) / total_height;
        float beta = (float)(y - t1.y) / segment_height; // be careful with divisions by zero 
        Vec2i A = t0 + (t2 - t0) * alpha;
        Vec2i B = t1 + (t2 - t1) * beta;
        line(A.x, y, B.x, y, image, red);
        //if (A.x > B.x) std::swap(A, B);
        //for (int j = A.x; j <= B.x; j++) {
        //    image.set(j, y, color); // attention, due to int casts t0.y+i != A.y 
        //}
    }
}

Vec3f barycentric(Vec3f* pts, Vec3f P) {
    Vec3f u = cross(Vec3f(pts[2][0] - pts[0][0], pts[1][0] - pts[0][0], pts[0][0] - P[0]), Vec3f(pts[2][1] - pts[0][1], pts[1][1] - pts[0][1], pts[0][1] - P[1]));
    /* `pts` and `P` has integer value as coordinates
       so `abs(u[2])` < 1 means `u[2]` is 0, that means
       triangle is degenerate, in this case return something with negative coordinates */
    if (std::abs(u[2]) < 1) return Vec3f(-1, 1, 1);
    return Vec3f(1.f - (u.x + u.y) / u.z, u.y / u.z, u.x / u.z);
}

void triangle(Model* model, Vec3f* pts, Vec2f* texture_coords, IShader& shader, TGAImage& image, TGAImage& zbuffer) {
    Vec2i bboxmin(image.get_width() - 1, image.get_height() - 1);
    Vec2i bboxmax(0, 0);
    Vec2i clamp(image.get_width() - 1, image.get_height() - 1);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
            bboxmin[j] = max(0, min(bboxmin[j], pts[i][j]));
            bboxmax[j] = min(clamp[j], max(bboxmax[j], pts[i][j]));
        }
    }
    Vec3f P;
    Vec2f Diffuse(0,0);
    for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++) {
        for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++) {

            Vec3f bc_screen = barycentric(pts, P);
            if (bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0) continue;
            P.z = 0;
            Diffuse[0] = 0;
            Diffuse[1] = 0;
            for (int i = 0; i < 3; i++) P.z += pts[i][2] * bc_screen[i];
            for (int i = 0; i < 3; i++) {
                Diffuse[0] += texture_coords[i][0] * bc_screen[i];
                Diffuse[1] += texture_coords[i][1] * bc_screen[i];
            }
            TGAColor color = model->diffuse(Diffuse);
            //TGAColor color = TGAColor(Diffuse[0], Diffuse[1], 0.5, 1);
            //float z = pts[0][2] * bc_screen.x + pts[1][2] * bc_screen.y + pts[2][2] * bc_screen.z;
            //float w = pts[0][3] * bc_screen.x + pts[1][3] * bc_screen.y + pts[2][3] * bc_screen.z;
            //float w = 1;
            //frag_depth = P.z;
            //float frag_depth = pts[2] * bc_clip;
            int frag_depth = int(max(0, min(255, P.z)));
            if (bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0 || zbuffer.get(P.x, P.y)[0] > frag_depth) continue;
            bool discard = shader.fragment(bc_screen, color);
            //if (zbuffer.get(P.x, P.y)[0] < frag_depth) {
            if (!discard) {
                zbuffer.set(P.x, P.y, TGAColor(frag_depth));
                image.set(P.x, P.y, color);
            }
        }
    }
}
