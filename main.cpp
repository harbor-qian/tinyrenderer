﻿#include <vector>
#include <windows.h>
#include "tgaimage.h"
#include "model.h"
#include "tinygl.h"

Model* model = NULL;
const int width  = 800;
const int height = 800;

Matrix ModelView;
Matrix Projection;
Matrix Viewport;

//Vec3f light_dir(0, 0, -1);
Vec3f light_dir(1, 1, 1);
//Vec3f       eye(1, 1, 3);
Vec3f       eye(0, 0, 3);
Vec3f    center(0, 0, 0);
Vec3f        up(0, 1, 0);

struct Shader : public IShader {
    mat<2, 3, float> varying_uv;  // same as above
    mat<4, 4, float> uniform_M;   //  Projection*ModelView
    mat<4, 4, float> uniform_MIT; // (Projection*ModelView).invert_transpose()

    virtual Vec4f vertex(int iface, int nthvert) {
        varying_uv.set_col(nthvert, model->uv(iface, nthvert));
        Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert)); // read the vertex from .obj file
        return Projection * ModelView * gl_Vertex; // transform it to screen coordinates
    }

    virtual bool fragment(Vec3f bar, TGAColor& color) {
        Vec2f uv = varying_uv * bar;
        Vec3f n = proj<3>(uniform_MIT * embed<4>(model->normal(uv))).normalize();
        Vec3f l = proj<3>(uniform_M * embed<4>(light_dir)).normalize();
        Vec3f r = (n * (n * l * 2.f) - l).normalize();   // reflected light
        float spec = pow(max(r.z, 0.0f), model->specular(uv));
        float diff = max(0.f, n * l);
        TGAColor c = model->diffuse(uv);
        color = c;
        for (int i = 0; i < 3; i++) color[i] = (float) min(5 + c[i] * (diff + 0.6 * spec), 255);
        return false;
    }
};
int main(int argc, char** argv) {
    if (2 > argc) {
        std::cerr << "Usage: " << argv[0] << " obj/model.obj" << std::endl;
        return 1;
    }

    TGAImage image(width, height, TGAImage::RGB);
    TGAImage zbuffer(width, height, TGAImage::GRAYSCALE);

    lookat(eye, center, up);
    viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4);
    //viewport(0, 0, width, height);
    projection(-1.f / (eye - center).norm());
    light_dir.normalize();
    //light_dir = proj<3>((Projection * ModelView * embed<4>(light_dir, 0.f))).normalize();
    for (int m = 1; m < argc; m++) {
        model = new Model(argv[m]);
        Shader shader;
        shader.uniform_M = Projection * ModelView;
        shader.uniform_MIT = (Projection * ModelView).invert_transpose();
        for (int i = 0; i < model->nfaces(); i++) {
            Vec4f clipc[3];
            Vec4f sc;
            int j = model->nverts();
            for (int j = 0; j < 3; j++) {
                Vec3f v = model->vert(i, j);
                clipc[j] = shader.vertex(i,j);
            }
            //Vec3f n = (world_coords[2] - world_coords[0]) ^ (world_coords[1] - world_coords[0]);
            //n.normalize();
            //float intensity = n * light_dir;
            //if (intensity > 0) {
            triangle(clipc, shader, image, zbuffer);
            //}
        }
    }
    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    zbuffer.flip_vertically();
    image.write_tga_file("output.tga");
    zbuffer.write_tga_file("zbuffer.tga");
    delete model;
    return 0;
}

