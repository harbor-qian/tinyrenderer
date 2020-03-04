#include <vector>
#include <windows.h>
#include "tgaimage.h"
#include "model.h"
#include "tinygl.h"

Model* model = NULL;
const int width  = 1024;
const int height = 1024;

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
    Vec3f varying_intensity; // written by vertex shader, read by fragment shader

    virtual Vec4f vertex(int iface, int nthvert) {
        varying_intensity[nthvert] = max(0.f, model->normal(iface, nthvert) * light_dir); // get diffuse lighting intensity
        Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert)); // read the vertex from .obj file
        return Viewport * Projection * ModelView * gl_Vertex; // transform it to screen coordinates
    }

    virtual bool fragment(Vec3f bar, TGAColor& color) {
        float intensity = varying_intensity * bar;   // interpolate intensity for the current pixel
        color = TGAColor(255, 255, 255) * intensity; // well duh
        return false;                              // no, we do not discard this pixel
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
    for (int m = 1; m < argc - 1; m++) {
        model = new Model(argv[m]);
        Shader shader;
        for (int i = 0; i < model->nfaces(); i++) {
            Vec3f screen_coords[3];
            //Vec3f world_coords[3];
            Vec2f uv[3];
            Vec4f sc;
            int j = model->nverts();
            for (int j = 0; j < 3; j++) {
                Vec3f v = model->vert(i, j);
                sc = Viewport * Projection * ModelView * Vec4f(v.x, v.y, v.z, 1);
                sc = shader.vertex(i,j);
                screen_coords[j] = Vec3f(sc.x / sc.w, sc.y / sc.w, sc.z / sc.w);
                //screen_coords[j] = Vec3f((v.x + 1.) * width / 2., (v.y + 1.) * height / 2., v.z);
                //world_coords[j] = v;
                uv[j] = model->uv(i, j);
            }
            //Vec3f n = (world_coords[2] - world_coords[0]) ^ (world_coords[1] - world_coords[0]);
            //n.normalize();
            //float intensity = n * light_dir;
            //if (intensity > 0) {
            triangle(model, screen_coords, uv, shader, image, zbuffer);
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

