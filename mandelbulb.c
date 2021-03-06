#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "mandelbulb.h"

struct Mandelbulb {
    int resolution;
    int (*get_iteration)(double, double, double);
    double A, B, C, D; // for quintic
};

struct VectorOfPoints {
    int size;
    int capacity;
    struct Point* data;
};

static int get_iteration_quintic(double x0, double y0, double z0) {
    double x = 0;
    double y = 0;
    double z = 0;
    double xn, yn, zn;
    int i;

    double A, B, C, D;
    A = B = C = D = 0;
    C = 2;

    for (i = 1; i<32; i=i+1) {
        double x2 = x*x, x3 = x*x2, x4 = x*x3;
        double y2 = y*y, y3 = y*y2, y4 = y*y3;
        double z2 = z*z, z3 = z*z2, z4 = z*z3;

        xn = x*x4 - 10*x3*(y2 + A*y*z + z2) +5*x*(y4 + B*y3*z + C*y2*z2 + B*y*z3 + z4)
            + D*x2*y*z*(y+z) + x0;
        yn = y*y4 - 10*y3*(z2 + A*x*z + x2) +5*y*(z4 + B*z3*x + C*z2*x2 + B*z*x3 + x4)
            + D*y2*z*x*(z+x) + y0;
        zn = z*z4 - 10*z3*(x2 + A*x*y + y2) +5*z*(x4 + B*x3*y + C*x2*y2 + B*x*y3 + y4)
            + D*z2*x*y*(x+y) + z0;

        if (sqrt(xn*xn+yn*yn+zn*zn) > 2) {
            return i;
        }
        x = xn; y = yn; z = zn;
    }
    return 0;
}

static int get_iteration_nine(double x0, double y0, double z0) {
    double x = 0;
    double y = 0;
    double z = 0;
    double xn, yn, zn;
    int i;

    for (i = 1; i<32; i=i+1) {
        double x2 = x*x, x3 = x*x2, x5 = x3*x2, x7 = x5*x2;
        double y2 = y*y, y3 = y*y2, y5 = y3*y2, y7 = y5*y2;
        double z2 = z*z, z3 = z*z2, z5 = z3*z2, z7 = z5*z2;

        xn = x7*x2 - 36*x7*(y2+z2) + 126*x5*(y2+z2)*(y2+z2)-84*x3*(y2+z2)*(y2+z2)*(y2+z2)
            +9*x*(y2+z2)*(y2+z2)*(y2+z2)*(y2+z2)+x0;
        yn = y7*y2 - 36*y7*(z2+x2) + 126*y5*(z2+x2)*(z2+x2)-84*y3*(z2+x2)*(z2+x2)*(z2+x2)
            +9*y*(z2+x2)*(z2+x2)*(z2+x2)*(z2+x2)+y0;
        zn = z7*z2 - 36*z7*(x2+y2) + 126*z5*(x2+y2)*(x2+y2)-84*z3*(x2+y2)*(x2+y2)*(x2+y2)
            +9*z*(x2+y2)*(x2+y2)*(x2+y2)*(x2+y2)+z0;

        if (sqrt(xn*xn+yn*yn+zn*zn) > 2) {
            return i;
        }
        x = xn; y = yn; z = zn;
    }
    return 0;
}

static int get_iteration_quadratic(double x0, double y0, double z0) {
    double x = 0;
    double y = 0;
    double z = 0;
    double xn, yn, zn;
    int i;

    for (i = 1; i<32; i=i+1) {
        xn = x*x-y*y-z*z+x0;
        yn = 2*x*z+y0;
        zn = 2*x*y+z0;

        if (sqrt(xn*xn+yn*yn+zn*zn) > 2) {
            return i;
        }
        x = xn; y = yn; z = zn;
    }
    return 0;
}

static int get_iteration_cubic(double x0, double y0, double z0) {
    double x = 0;
    double y = 0;
    double z = 0;
    double xn, yn, zn;
    int i;

    for (i = 1; i<32; i=i+1) {
        xn = x*x*x-3*x*(y*y+z*z)+x0;
        yn = -y*y*y+3*y*x*x-y*z*z+y0;
        zn = z*z*z-3*z*x*x+z*y*y+z0;

        if (sqrt(xn*xn+yn*yn+zn*zn) > 2) {
            return i;
        }
        x = xn; y = yn; z = zn;
    }
    return 0;
}

struct Mandelbulb* mandelbulb_quintic_new(int resolution, double A, double B, double C, double D) {
    struct Mandelbulb* m = calloc(1, sizeof(struct Mandelbulb));
    m->resolution = resolution;
    m->get_iteration = get_iteration_quintic;
    m->A = A;
    m->B = B;
    m->C = C;
    m->D = D;
    return m;
}

struct Mandelbulb* mandelbulb_nine_new(int resolution) {
    struct Mandelbulb* m = calloc(1, sizeof(struct Mandelbulb));
    m->resolution = resolution;
    m->get_iteration = get_iteration_nine;
    return m;
}

struct Mandelbulb* mandelbulb_quadratic_new(int resolution) {
    struct Mandelbulb* m = calloc(1, sizeof(struct Mandelbulb));
    m->resolution = resolution;
    m->get_iteration = get_iteration_quadratic;
    return m;
}

struct Mandelbulb* mandelbulb_cubic_new(int resolution) {
    struct Mandelbulb* m = calloc(1, sizeof(struct Mandelbulb));
    m->resolution = resolution;
    m->get_iteration = get_iteration_cubic;
    return m;
}

void mandelbulb_free(struct Mandelbulb* mandelbulb) {
    free(mandelbulb);
}

struct Point* mandelbulb_surface_points(struct Mandelbulb* mandelbulb, long* size) {
    int stride = mandelbulb->resolution;
    long width = stride;
    long height = stride;
    long depth = stride;
    unsigned long long* mask = calloc(1, width*height*depth/8);
#pragma omp parallel for
    for (long zbyte = 0; zbyte < depth/64; zbyte=zbyte+1) {
        for (long z = zbyte*64; z < (zbyte+1)*64; z=z+1) {
     //    for (int z = 0; z < depth; z=z+1) {
            for (long y = 0; y < height; y=y+1) {
                for (long x = 0; x < width; x=x+1) {
                    double v[3] = {x,y,z};
                    for (int i = 0; i < 3; i=i+1) {
                        v[i] /= stride;
                        v[i] -= 0.5;
                        v[i] *= 4;
                    }

                    int it = mandelbulb->get_iteration(v[0], v[1], v[2]);
                    if (it == 0) {
                        long byte = (z*width*height+y*width+x)/64;
                        long bit = (z*width*height+y*width+x)%64;
                        mask[byte] |= (1ULL << bit);
                    }
                }
            }
        }
    }
    struct VectorOfPoints Z[2048];
    memset(Z, 0, sizeof(Z));
#pragma omp parallel for
    for (long z = 1; z < depth-1; z=z+1) {
        for (long y = 1; y < height-1; y=y+1) {
            for (long x = 1; x < width-1; x=x+1) {
                int byte = (z*width*height+y*width+x)/64;
                int bit = (z*width*height+y*width+x)%64;
                if (! (mask[byte] & (1ULL << bit))) {
                    continue;
                }
                for (long i = -1; i <= 1; i=i+1) {
                    for (long j = -1; j <= 1; j=j+1) {
                        for (long k = -1; k <= 1; k=k+1) {
                            long byte = ((z+i)*width*height+(y+j)*width+(x+k))/64;
                            long bit = ((z+i)*width*height+(y+j)*width+(x+k))%64;
                            if (! (mask[byte] & (1ULL << bit))) {
                                if (Z[z].size >= Z[z].capacity) {
                                    Z[z].capacity = (Z[z].capacity+1)*2;
                                    Z[z].data = realloc(Z[z].data, Z[z].capacity*sizeof(struct Point));
                                }
                                Z[z].data[Z[z].size].v[0] = (x/(double)width)-0.5;
                                Z[z].data[Z[z].size].v[1] = (y/(double)height)-0.5;
                                Z[z].data[Z[z].size].v[2] = (z/(double)depth)-0.5;
                                Z[z].size++;
                                goto done;
                            }
                        }
                    }
                }
done:
            ;
            }
        }
    }
    fprintf(stderr, "builing surface\n");
    int boundary= 0;
    for (int z = 1; z < depth-1; z=z+1) {
        boundary += Z[z].size;
    }
    struct Point* points = malloc(boundary * sizeof(struct Point));
    int npoints = 0;
    for (int z = 1; z < depth-1; z=z+1) {
        memcpy(&points[npoints], Z[z].data, Z[z].size*sizeof(struct Point));
        npoints += Z[z].size;
        free(Z[z].data);
    }
    fprintf(stderr, "boundary: %d\n", boundary);
    free(mask);

    *size = boundary;

    return points;
}

struct VertexInfo* mandelbulb_surface_vertices(struct Mandelbulb* mandelbulb, long* size) {
    long npoints;
    struct Point* points = mandelbulb_surface_points(mandelbulb, &npoints);
    float d = 1.5/mandelbulb->resolution;
    long nvertex = npoints*6*2*3;
    struct VertexInfo* vertex_data = malloc(nvertex*sizeof(struct VertexInfo));

    int j = 0;
    for (int i = 0; i < npoints; i=i+1) {
        float x0 = points[i].v[0], y0 = points[i].v[1], z0 = points[i].v[2];
        x0 *= 2; y0 *= 2; z0 *= 2;
        struct VertexInfo vd[] = {
        // up
        {{x0-d, y0-d, z0+d}, {1.0f, 1.0f, 0.0f}},
        {{x0-d, y0+d, z0+d}, {1.0f, 1.0f, 0.0f}},
        {{x0+d, y0-d, z0+d}, {1.0f, 1.0f, 0.0f}},

        {{x0+d, y0+d, z0+d}, {1.0f, 1.0f, 0.0f}},
        {{x0-d, y0+d, z0+d}, {1.0f, 1.0f, 0.0f}},
        {{x0+d, y0-d, z0+d}, {1.0f, 1.0f, 0.0f}},

        // down
        {{x0-d, y0-d, z0-d}, {1.0f, 1.0f, 0.0f}},
        {{x0-d, y0+d, z0-d}, {1.0f, 1.0f, 0.0f}},
        {{x0+d, y0-d, z0-d}, {1.0f, 1.0f, 0.0f}},

        {{x0+d, y0+d, z0-d}, {1.0f, 1.0f, 0.0f}},
        {{x0-d, y0+d, z0-d}, {1.0f, 1.0f, 0.0f}},
        {{x0+d, y0-d, z0-d}, {1.0f, 1.0f, 0.0f}},

        // left
        {{x0-d, y0-d, z0-d}, {1.0f, 1.0f, 0.0f}},
        {{x0-d, y0+d, z0-d}, {1.0f, 1.0f, 0.0f}},
        {{x0-d, y0-d, z0+d}, {1.0f, 1.0f, 0.0f}},

        {{x0-d, y0+d, z0-d}, {1.0f, 1.0f, 0.0f}},
        {{x0-d, y0+d, z0+d}, {1.0f, 1.0f, 0.0f}},
        {{x0-d, y0-d, z0+d}, {1.0f, 1.0f, 0.0f}},

        // right
        {{x0+d, y0-d, z0-d}, {1.0f, 1.0f, 0.0f}},
        {{x0+d, y0+d, z0-d}, {1.0f, 1.0f, 0.0f}},
        {{x0+d, y0-d, z0+d}, {1.0f, 1.0f, 0.0f}},

        {{x0+d, y0+d, z0-d}, {1.0f, 1.0f, 0.0f}},
        {{x0+d, y0+d, z0+d}, {1.0f, 1.0f, 0.0f}},
        {{x0+d, y0-d, z0+d}, {1.0f, 1.0f, 0.0f}},

        // front
        {{x0-d, y0-d, z0-d}, {1.0f, 1.0f, 0.0f}},
        {{x0+d, y0-d, z0-d}, {1.0f, 1.0f, 0.0f}},
        {{x0-d, y0-d, z0+d}, {1.0f, 1.0f, 0.0f}},

        {{x0+d, y0-d, z0-d}, {1.0f, 1.0f, 0.0f}},
        {{x0+d, y0-d, z0+d}, {1.0f, 1.0f, 0.0f}},
        {{x0-d, y0-d, z0+d}, {1.0f, 1.0f, 0.0f}},

        // rear
        {{x0-d, y0+d, z0-d}, {1.0f, 1.0f, 0.0f}},
        {{x0+d, y0+d, z0-d}, {1.0f, 1.0f, 0.0f}},
        {{x0-d, y0+d, z0+d}, {1.0f, 1.0f, 0.0f}},

        {{x0+d, y0+d, z0-d}, {1.0f, 1.0f, 0.0f}},
        {{x0+d, y0+d, z0+d}, {1.0f, 1.0f, 0.0f}},
        {{x0-d, y0+d, z0+d}, {1.0f, 1.0f, 0.0f}}
        };

        memcpy(&vertex_data[j], vd,6*2*3*sizeof(struct VertexInfo));
        j+=6*2*3;
    }

    *size = nvertex;
    free(points);

    return vertex_data;
}
