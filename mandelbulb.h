#pragma once

struct Mandelbulb;

struct Point {
    float v[3];
};

struct VertexInfo {
    float position[3];
    float color[3];
};

struct Mandelbulb* mandelbulb_quintic_new(int resolution, double A, double B, double C, double D);
struct Mandelbulb* mandelbulb_nine_new(int resolution);
struct Mandelbulb* mandelbulb_quadratic_new(int resolution);
struct Mandelbulb* mandelbulb_cubic_new(int resolution);

void mandelbulb_free(struct Mandelbulb* mandelbulb);

struct Point* mandelbulb_surface_points(struct Mandelbulb* mandelbulb, long* size);
struct VertexInfo* mandelbulb_surface_vertices(struct Mandelbulb* mandelbulb, long* size);
