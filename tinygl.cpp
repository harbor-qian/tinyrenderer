#include "tinygl.h"
#include <cmath>
//FILE* fp = fopen("z.txt", "w");


IShader::~IShader() {}

void viewport(int x, int y, int w, int h) {
    Viewport = Matrix::identity();
    Viewport[0][0] = w / 2.f;
    Viewport[0][3] = x + w / 2.f;
    Viewport[1][1] = h / 2.f;
    Viewport[1][3] = y + h / 2.f;
    //Viewport[2][2] = -255/2.f;
    //Viewport[2][3] = 255/2.f;
    // hack, doesn't use this 
    float n = 1;
    float f = 3;
    Viewport[2][2] = (f-n)/2.f;
    Viewport[2][3] = (f+n)/2.f;
}

void projection(float fovy, float aspect, float n, float f) {  // n, f>0
    Projection = Matrix::identity();
    //Projection[3][2] = coeff;
#ifdef ORTHO_PROJ
    Projection[0][0] = fovy / aspect;
    Projection[1][1] = fovy;
    Projection[2][2] = -2/(f-n);
    Projection[2][3] = -(f + n) / (f - n);
#else
    Projection[0][0] = fovy/aspect;  // tan(fovy/2)*n = t;
    Projection[1][1] = fovy;
    Projection[2][2] = -(f+n)/(f-n);
    Projection[2][3] = -2*f*n/(f-n);
    Projection[3][2] = -1;
    Projection[3][3] = 0;
#endif
}

void lookat(Vec3f eye, Vec3f center, Vec3f up) {
    Vec3f z = (eye - center).normalize();  // -forward
    Vec3f x = cross(up, z).normalize();    // side
    Vec3f y = cross(z, x).normalize();     // up
    ModelView = Matrix::identity();
    Matrix Translate = Matrix::identity();
    Vec3f tmp;
    tmp[0] = (eye * x);
    tmp[1] = (eye * y);
    tmp[2] = (eye * z);
    eye[2] += 2;  //hack the obj position as the origin is at (0,0,0)
    for (int i = 0; i < 3; i++) {
        ModelView[0][i] = x[i];
        ModelView[1][i] = y[i];
        ModelView[2][i] = z[i];
    //    ModelView[i][3] = -tmp[i];
        Translate[i][3] = -eye[i];
    }
    // (Translate * ModelView)-1 = (ModelView)T * Translate(-e)
    ModelView = ModelView * Translate;
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

float EdgeFunc(const Vec3f& p0, const Vec3f& p1, const Vec3f& p2) {
    return ((p2.x - p0.x) * (p1.y - p0.y) - (p2.y - p0.y) * (p1.x - p0.x));
} // note that the result of edge function could be represent as area as well.


//void triangle(Vec4f* clipc, IShader& shader, TGAImage& image, TGAImage& zbuffer) {
void triangle(Vec4f* clipc, IShader& shader, TGAImage& image, float* zbuffer) {
    Vec2i bboxmin(image.get_width() - 1, image.get_height() - 1);
    Vec2i bboxmax(0, 0);
    Vec2i clamp(image.get_width() - 1, image.get_height() - 1);
    Vec3f pts[3];
    //fprintf(fp, "new triangle\n");

    for (int i = 0; i < 3; i++) {
        pts[i] = proj<3>(Viewport * clipc[i] / clipc[i][3]);
        pts[i][2] = -clipc[i][3];
    }
    float area = EdgeFunc(pts[0], pts[1], pts[2]);
    Vec3f P10 = pts[1]-pts[0];
    Vec3f P20 = pts[2]-pts[0];
    // z component of (P10 x P20), ccw  |P10.x, P10.y|, same as edge function
    //                                  |P20.x, P20.y|
    // = area, ccw
    bool frontface = P10.x * P20.y - P10.y * P20.x > 0; //>0, P1 on P20's right, frontface; <0, on P20's left
    //if (!frontface) return;
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
            bboxmin[j] = max(0, min(bboxmin[j], pts[i][j]));
            bboxmax[j] = min(clamp[j], max(bboxmax[j], pts[i][j]));
        }
    }
    Vec3i P;
    for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++) {
        for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++) {

            Vec3f bc_screen = barycentric(pts, P);
            //Vec3f bc_screen;
            //bc_screen.x= EdgeFunc(pts[1], pts[2], P)/area;
            //bc_screen.y= EdgeFunc(pts[2], pts[0], P)/area;
            //bc_screen.z= EdgeFunc(pts[0], pts[1], P)/area;
            // u/z, v/z interp, calculate 1/z barycentric first
            Vec3f bc_clip = Vec3f(bc_screen.x / clipc[0][3], bc_screen.y / clipc[1][3], bc_screen.z / clipc[2][3]);
            float w = 1/(bc_clip.x + bc_clip.y + bc_clip.z);
            //fprintf(fp, "%f\n", w);

            float frag_depth = w;
            bc_clip = bc_clip*w;// two "-" cancel out
            
#ifdef ORTHO_PROJ
            frag_depth = 0;
            for (int i = 0; i < 3; i++) frag_depth += clipc[i][2] * bc_clip[i];
#endif
            //float n = 2;
            //float f = 5;
            //P.z = (f-P.z) / (f - n);

            //int frag_depth = int(max(0, min(255, P.z*255.0f)));
            //if (bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0 || zbuffer.get(P.x, P.y)[0] > frag_depth) continue;
            //printf("%f\n", frag_depth);
            if (bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0 || zbuffer[P.x + P.y * image.get_width()] < frag_depth) continue;
            TGAColor color;
            bool discard = shader.fragment(bc_clip, color);
            //bool discard = shader.fragment(bc_screen, color);
            if (!discard) {
                zbuffer[P.x + P.y * image.get_width()] = frag_depth;
                //zbuffer.set(P.x, P.y, TGAColor(frag_depth));
                image.set(P.x, P.y, color);
            }
        }
    }
}
