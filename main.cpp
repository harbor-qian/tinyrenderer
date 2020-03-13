#include <vector>
#include <windows.h>
#include <climits>
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
//Vec3f       eye(1, 1, 0);
Vec3f       eye(0, 0, 0);
//Vec3f       eye(0, 0, 3);
Vec3f    center(0, 0, -10);
//Vec3f    center(1, 1, -3);
Vec3f        up(0, 1, 0);

struct Shader : public IShader {
    mat<2, 3, float> varying_uv;  // same as above
    mat<3, 3, float> varying_nrm;
    mat<3, 3, float> ndc_tri;
    mat<4, 4, float> uniform_M;   //  Projection*ModelView
    mat<4, 4, float> uniform_MIT; // (Projection*ModelView).invert_transpose()


    virtual Vec4f vertex(int iface, int nthvert) {
        varying_nrm.set_col(nthvert, proj<3>((Projection * ModelView).invert_transpose() * embed<4>(model->normal(iface, nthvert), 0.f)));

        Vec4f gl_Vertex = ModelView *embed<4>(model->vert(iface, nthvert)); // read the vertex from .obj file
        gl_Vertex = Projection * ModelView *embed<4>(model->vert(iface, nthvert)); // read the vertex from .obj file
        varying_uv.set_col(nthvert, model->uv(iface, nthvert));
        //varying_uv.set_col(nthvert, model->uv(iface, nthvert)/gl_Vertex[3]);
        ndc_tri.set_col(nthvert, proj<3>(gl_Vertex / gl_Vertex[3]));
        return gl_Vertex; // transform it to screen coordinates
    }

    virtual bool fragment(Vec3f bar, TGAColor& color) {
        Vec3f bn = (varying_nrm * bar).normalize();
        Vec2f uv = varying_uv * bar;

        float diff = max(0.f, bn * light_dir);
        mat<3, 3, float> A;
        A[0] = ndc_tri.col(1) - ndc_tri.col(0);
        A[1] = ndc_tri.col(2) - ndc_tri.col(0);
        A[2] = bn;

        mat<3, 3, float> AI = A.invert();

        Vec3f i = AI * Vec3f(varying_uv[0][1] - varying_uv[0][0], varying_uv[0][2] - varying_uv[0][0], 0);
        Vec3f j = AI * Vec3f(varying_uv[1][1] - varying_uv[1][0], varying_uv[1][2] - varying_uv[1][0], 0);

        mat<3, 3, float> B;
        B.set_col(0, i.normalize());
        B.set_col(1, j.normalize());
        B.set_col(2, bn);

        Vec3f n = (B * model->normal(uv)).normalize();

        color = model->diffuse(uv) * diff;

        //Vec3f n = proj<3>(uniform_MIT * embed<4>(model->normal(uv))).normalize();
        Vec3f l = proj<3>(uniform_M * embed<4>(light_dir)).normalize();
        Vec3f r = (n * (n * l * 2.f) - l).normalize();   // reflected light
        //float diff = max(0.f, n * light_dir);
        float spec = pow(max(r.z, 0.0f), model->specular(uv));
        diff = max(0.f, n * l);
        color = model->diffuse(uv);
        //for (int i = 0; i < 3; i++) color[i] = (float) min(5 + color[i] * (diff + 0.6 * spec), 255);
        return false;
    }
};
int main(int argc, char** argv) {
    if (2 > argc) {
        std::cerr << "Usage: " << argv[0] << " obj/model.obj" << std::endl;
        return 1;
    }

    TGAImage image(width, height, TGAImage::RGB);
    //TGAImage zbuffer(width, height, TGAImage::GRAYSCALE);
    float* zbuffer = new float[width * height];
    for (int i = 0; i < width * height; i++)
        zbuffer[i] = FLT_MAX;

    lookat(eye, center, up);
    //float n = eye.z - 1;
    float n = 1;
    float f = 3;
    float t = 1;   //max yG
    float r = 1;   //max x
    float aspect = r/t;
    float fovy = n / t;
    projection(fovy,aspect, n, f);
    //projection(-1.f / (eye - center).norm());
    //viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4);
    viewport(0, 0, width, height);
    light_dir.normalize();
    //light_dir = proj<3>((Projection * ModelView * embed<4>(light_dir, 0.f))).normalize();
    for (int m = 1; m < argc; m++) {
        model = new Model(argv[m]);
        Shader shader;
        shader.uniform_M = Projection * ModelView;
        shader.uniform_MIT = (Projection * ModelView).invert_transpose();
        for (int i = 0; i < model->nfaces(); i++) {
            Vec4f clipc[3];
            Vec3f clip[3];
            Vec4f sc;
            int j = model->nverts();
            for (int j = 0; j < 3; j++) {
                Vec3f v = model->vert(i, j);
                clipc[j] = shader.vertex(i, j);
                //clip[j] = proj<3>(clipc[j]);
                //clip[j] = proj<3>(ModelView * embed<4>(v));
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
    //zbuffer.flip_vertically();
    image.write_tga_file("output.tga");
    //zbuffer.write_tga_file("zbuffer.tga");
    delete model;
    return 0;
}

