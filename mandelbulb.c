#include <gtk/gtk.h>
#include <math.h>
#include <limits.h>
#include <stddef.h>

#include <epoxy/gl.h>

struct Point {
    float v[3];
};

struct VertexInfo {
    float position[3];
    float color[3];
};

struct VectorOfPoints {
    int size;
    int capacity;
    struct Point* data;
};

struct App {
    GtkWidget* window;
    GtkWidget* glarea;
    double tx, ty, tz;

    guint program;
    guint mvp_location;
    guint position_location;
    guint color_location;
    guint vao;

    float mvp[16];

    struct VertexInfo* vertex_data;
    int nvertex;
};

static void
compute_mvp (float *res,
             float  phi,
             float  theta,
             float  psi)
{ 
    float x = phi * (G_PI / 180.f);
    float y = theta * (G_PI / 180.f);
    float z = psi * (G_PI / 180.f);
    float c1 = cosf (x), s1 = sinf (x);
    float c2 = cosf (y), s2 = sinf (y);
    float c3 = cosf (z), s3 = sinf (z);
    float c3c2 = c3 * c2;
    float s3c1 = s3 * c1;
    float c3s2s1 = c3 * s2 * s1;
    float s3s1 = s3 * s1;
    float c3s2c1 = c3 * s2 * c1;
    float s3c2 = s3 * c2;
    float c3c1 = c3 * c1;
    float s3s2s1 = s3 * s2 * s1;
    float c3s1 = c3 * s1;
    float s3s2c1 = s3 * s2 * c1;
    float c2s1 = c2 * s1;
    float c2c1 = c2 * c1;

    /* apply all three Euler angles rotations using the three matrices:
    *
    * ⎡  c3 s3 0 ⎤ ⎡ c2  0 -s2 ⎤ ⎡ 1   0  0 ⎤
    * ⎢ -s3 c3 0 ⎥ ⎢  0  1   0 ⎥ ⎢ 0  c1 s1 ⎥
    * ⎣   0  0 1 ⎦ ⎣ s2  0  c2 ⎦ ⎣ 0 -s1 c1 ⎦
    */
    res[0] = c3c2;  res[4] = s3c1 + c3s2s1;  res[8] = s3s1 - c3s2c1; res[12] = 0.f;
    res[1] = -s3c2; res[5] = c3c1 - s3s2s1;  res[9] = c3s1 + s3s2c1; res[13] = 0.f;
    res[2] = s2;    res[6] = -c2s1;         res[10] = c2c1;          res[14] = 0.f;
    res[3] = 0.f;   res[7] = 0.f;           res[11] = 0.f;           res[15] = 1.f;
}

static gboolean key_press_event(GtkWidget* widget, GdkEventKey* event, struct App* app) 
{
    if (event->type == GDK_KEY_PRESS) {
        if (event->keyval == GDK_KEY_Down) {
            app->tx -= 2*(2*M_PI)/36.0;
            compute_mvp(app->mvp, app->tx, app->ty, app->tz);
            gtk_widget_queue_draw(app->glarea);
        } else if (event->keyval == GDK_KEY_Up) {
            app->tx += 2*(2*M_PI)/36.0;
            compute_mvp(app->mvp, app->tx, app->ty, app->tz);
            gtk_widget_queue_draw(app->glarea);
        } else if (event->keyval == GDK_KEY_Left) {
            app->ty += 2*(2*M_PI)/36.0;
            compute_mvp(app->mvp, app->tx, app->ty, app->tz);
            gtk_widget_queue_draw(app->glarea);
        } else if (event->keyval == GDK_KEY_Right) {
            app->ty -= 2*(2*M_PI)/36.0;
            compute_mvp(app->mvp, app->tx, app->ty, app->tz);
            gtk_widget_queue_draw(app->glarea);
        } else {
            printf("%x\n", event->keyval);
        }
    }
    return TRUE;
}

int get_iteration(double x0, double y0, double z0) {
//    int n = 8;
    double x = 0;
    double y = 0;
    double z = 0;
    double xn, yn, zn;
    int i;

    double A, B, C, D;
    A = B = C = D = 0;
    C = 2;

    for (i = 1; i<32; i=i+1) {
        // wikipedia

        /*
        double r = sqrt(x*x+y*y+z*z);
        double theta = atan(y/x);
        double phi = atan(z/sqrt(x*x+y*y));

        double rn = r*r*r*r*r*r*r*r;
        xn = rn*cos(n*theta)*cos(n*phi) + x0;
        yn = rn*sin(n*theta)*cos(n*phi) + y0;
        zn = rn*sin(n*phi) + z0;
        */

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

static void calc_surface(struct App* app) {
    int stride = 1024; // 128; // 1024;
    int width = stride;
    int height = stride;
    int depth = stride;
    unsigned long long* mask = calloc(1, width*height*depth);
#pragma omp parallel for
    for (int zbyte = 0; zbyte < depth/64; zbyte=zbyte+1) {
        for (int z = zbyte*64; z < (zbyte+1)*64; z=z+1) {
     //    for (int z = 0; z < depth; z=z+1) {
            for (int y = 0; y < height; y=y+1) {
                for (int x = 0; x < width; x=x+1) {
                    double v[3] = {x,y,z};
                    for (int i = 0; i < 3; i=i+1) {
                        v[i] /= stride;
                        v[i] -= 0.5;
                        v[i] *= 4;
                    }

                    int it = get_iteration(v[0], v[1], v[2]);
                    if (it == 0) {
                        int byte = (z*width*height+y*width+x)/64;
                        int bit = (z*width*height+y*width+x)%64;
                        mask[byte] |= (1ULL << bit);
                    }
                }
            }
        }
    }

    struct VectorOfPoints Z[1024];
    memset(Z, 0, sizeof(Z));
#pragma omp parallel for
    for (int z = 1; z < depth-1; z=z+1) {
        for (int y = 1; y < height-1; y=y+1) {
            for (int x = 1; x < width-1; x=x+1) {
                int byte = (z*width*height+y*width+x)/64;
                int bit = (z*width*height+y*width+x)%64;
                if (! (mask[byte] & (1ULL << bit))) {
                    continue;
                }
                for (int i = -1; i <= 1; i=i+1) {
                    for (int j = -1; j <= 1; j=j+1) {
                        for (int k = -1; k <= 1; k=k+1) {
                            int byte = ((z+i)*width*height+(y+j)*width+(x+k))/64;
                            int bit = ((z+i)*width*height+(y+j)*width+(x+k))%64;
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
    }
    fprintf(stderr, "boundary: %d\n", boundary);
    free(mask);

    float d = 1.5/stride;
    int nvertex = npoints*6*2*3;
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
        {{x0-d, y0-d, z0+d}, {1.0f, 1.0f, 0.0f}},
        {{x0-d, y0+d, z0+d}, {1.0f, 1.0f, 0.0f}},
        {{x0+d, y0-d, z0+d}, {1.0f, 1.0f, 0.0f}},

        {{x0+d, y0+d, z0+d}, {1.0f, 1.0f, 0.0f}},
        {{x0-d, y0+d, z0+d}, {1.0f, 1.0f, 0.0f}},
        {{x0+d, y0-d, z0+d}, {1.0f, 1.0f, 0.0f}},

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
        {{x0-d, y0-d, z0-d}, {1.0f, 1.0f, 0.0f}},
        {{x0+d, y0-d, z0-d}, {1.0f, 1.0f, 0.0f}},
        {{x0-d, y0-d, z0+d}, {1.0f, 1.0f, 0.0f}},

        {{x0+d, y0-d, z0-d}, {1.0f, 1.0f, 0.0f}},
        {{x0+d, y0-d, z0+d}, {1.0f, 1.0f, 0.0f}},
        {{x0-d, y0-d, z0+d}, {1.0f, 1.0f, 0.0f}}
        };

        memcpy(&vertex_data[j], vd,6*2*3*sizeof(struct VertexInfo));
        j+=6*2*3;
    }

    app->vertex_data = vertex_data;
    app->nvertex = nvertex;
    free(points);
}

static gboolean
gl_render (GtkGLArea *area, GdkGLContext *context, struct App* app)
{   
    /* clear the viewport; the viewport is automatically resized when
     * the GtkGLArea gets a new size allocation
     */
    glClearColor (0.5, 0.5, 0.5, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    /* draw our object */
    if (app->program == 0 || app->vao == 0) {
        return FALSE;
    }

    /* load our program */
    glUseProgram (app->program);

    /* update the "mvp" matrix we use in the shader */
    glUniformMatrix4fv (app->mvp_location, 1, GL_FALSE, &(app->mvp[0]));

    /* use the buffers in the VAO */
    glBindVertexArray (app->vao);

    /* draw the three vertices as a triangle */
    glDrawArrays (GL_TRIANGLES, 0, app->nvertex);

    /* we finished using the buffers and program */
    glBindVertexArray (0);
    glUseProgram (0);

    /* flush the contents of the pipeline */
    glFlush ();

    return FALSE;
}

static guint
create_shader (int          shader_type,
               const char  *source,
               guint       *shader_out)
{ 
    guint shader = glCreateShader (shader_type);
    glShaderSource (shader, 1, &source, NULL);
    glCompileShader (shader);

    int status;
    glGetShaderiv (shader, GL_COMPILE_STATUS, &status);
    if (status == GL_FALSE)
    { 
        int log_len;
        glGetShaderiv (shader, GL_INFO_LOG_LENGTH, &log_len);
        
        char *buffer = g_malloc (log_len + 1);
        glGetShaderInfoLog (shader, log_len, NULL, buffer);
        
        fprintf(stderr,
                "Compilation failure in %s shader: %s",
                shader_type == GL_VERTEX_SHADER ? "vertex" : "fragment",
                buffer);
        
        g_free (buffer);

        glDeleteShader (shader);
        shader = 0;
    }

    if (shader_out != NULL) {
        *shader_out = shader;
    }

    return shader != 0;
}

static gboolean
init_shaders (struct App* app)
{
    const char* vertex_shader = " \n \
#version 150 \n \
\n \
in vec3 position; \n \
in vec3 color; \n \
\n \
uniform mat4 mvp; \n \
\n \
smooth out vec4 vertexColor; \n \
\n \
void main() { \n \
  gl_Position = mvp * vec4(position, 1.0); \n \
\n \
  vec3 linear = pow(color, vec3(1.0 / 2.2)); \n \
  float dst = distance(vec4(0, 0, 2, 1), gl_Position) * 0.4; \n \
  vec3 scaled = linear * vec3(dst); \n \
\n \
  vertexColor = vec4(pow(scaled, vec3(2.2)), 1.0); \n \
} \n \
    ";
    const char* fragment_shader = "\n \
#version 150 \n \
smooth in vec4 vertexColor; \n \
\n \
out vec4 outputColor; \n \
\n \
void main() { \n \
  outputColor = vertexColor; \n \
}";
    guint program = 0;
    guint mvp_location = 0;
    guint vertex = 0, fragment = 0;
    guint position_location = 0;
    guint color_location = 0;

    /* load the vertex shader */
    create_shader (GL_VERTEX_SHADER, vertex_shader, &vertex);
    if (vertex == 0) {
        goto out;
    }

    /* load the fragment shader */
    create_shader (GL_FRAGMENT_SHADER, fragment_shader, &fragment);
    if (fragment == 0) {
        goto out;
    }

    /* link the vertex and fragment shaders together */
    program = glCreateProgram ();
    glAttachShader (program, vertex);
    glAttachShader (program, fragment);
    glLinkProgram (program);

    int status = 0;
    glGetProgramiv (program, GL_LINK_STATUS, &status);
    if (status == GL_FALSE)
    {
        int log_len = 0;
        glGetProgramiv (program, GL_INFO_LOG_LENGTH, &log_len);

        char *buffer = g_malloc (log_len + 1);
        glGetProgramInfoLog (program, log_len, NULL, buffer);

        fprintf(stderr, "Linking failure in program: %s", buffer);

        g_free (buffer);

        glDeleteProgram (program);
        program = 0;

        goto out;
    }

    /* get the location of the "mvp" uniform */
    mvp_location = glGetUniformLocation (program, "mvp");

    /* get the location of the "position" and "color" attributes */
    position_location = glGetAttribLocation (program, "position");
    color_location = glGetAttribLocation (program, "color");

    /* the individual shaders can be detached and destroyed */
    glDetachShader (program, vertex);
    glDetachShader (program, fragment);

out:
    if (vertex != 0) {
        glDeleteShader (vertex);
    }
    if (fragment != 0) {
        glDeleteShader (fragment);
    }

    app->program = program;
    app->mvp_location = mvp_location;
    app->position_location = position_location;
    app->color_location = color_location;

    return program != 0;
}

static void init_buffers(struct App* app) {
    fprintf(stderr, "calc surface\n");
    calc_surface(app);
    fprintf(stderr, "calc surface done\n");
    fprintf(stderr, "nvertex %d\n", app->nvertex);

    guint vao, buffer;
    guint color_index = app->color_location;
    guint position_index = app->position_location;
    struct VertexInfo* vertex_data = app->vertex_data;
    int nvertex = app->nvertex;

    /* we need to create a VAO to store the other buffers */
    glGenVertexArrays (1, &vao);
    glBindVertexArray (vao);

    /* this is the VBO that holds the vertex data */
    glGenBuffers (1, &buffer);
    glBindBuffer (GL_ARRAY_BUFFER, buffer);
    glBufferData (GL_ARRAY_BUFFER, nvertex*sizeof(struct VertexInfo), vertex_data, GL_STATIC_DRAW);

    /* enable and set the position attribute */
    glEnableVertexAttribArray (position_index);
    glVertexAttribPointer (
        position_index, 3, GL_FLOAT, GL_FALSE,
        sizeof (struct VertexInfo),
        (GLvoid *) (G_STRUCT_OFFSET (struct VertexInfo, position)));

    /* enable and set the color attribute */
    glEnableVertexAttribArray (color_index);
    glVertexAttribPointer (
        color_index, 3, GL_FLOAT, GL_FALSE,
        sizeof (struct VertexInfo),
        (GLvoid *) (G_STRUCT_OFFSET (struct VertexInfo, color)));

    /* reset the state; we will re-enable the VAO when needed */
    glBindBuffer (GL_ARRAY_BUFFER, 0);
    glBindVertexArray (0);

    /* the VBO is referenced by the VAO */
    glDeleteBuffers (1, &buffer);

    app->vao = vao;
}

static void gl_init(GtkWidget* widget, struct App* app) {
    gtk_gl_area_make_current (GTK_GL_AREA(widget));

    if (gtk_gl_area_get_error (GTK_GL_AREA(widget)) != NULL) {
        return;
    }

    if (!init_shaders(app))
    {
        fprintf(stderr, "cannot init shaders\n");
        return;
    }

    init_buffers(app);
}

static void
init_mvp (float *res)
{
    /* initialize a matrix as an identity matrix */
    res[0] = 1.f; res[4] = 0.f;  res[8] = 0.f; res[12] = 0.f;
    res[1] = 0.f; res[5] = 1.f;  res[9] = 0.f; res[13] = 0.f;
    res[2] = 0.f; res[6] = 0.f; res[10] = 1.f; res[14] = 0.f;
    res[3] = 0.f; res[7] = 0.f; res[11] = 0.f; res[15] = 1.f;
}

static void close_window(GtkWidget* window, struct App* app)
{
    /* TODO: cleanup */
    gtk_main_quit();
}

int main(int argc, char** argv) {
    struct App app;
    memset(&app, 0, sizeof(struct App));
    GtkWidget* window;
    GtkWidget* glarea;
    gtk_init(&argc, &argv);

    init_mvp (app.mvp);

    window = gtk_window_new (GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(window), "Mandelbulb");
    gtk_window_set_default_size(GTK_WINDOW(window), 800, 800);
    glarea = gtk_gl_area_new();
    gtk_gl_area_set_has_depth_buffer(GTK_GL_AREA(glarea), TRUE);

    gtk_container_add(GTK_CONTAINER(window), glarea);

    app.window = window;
    app.glarea = glarea;

    g_signal_connect(window, "key-press-event", G_CALLBACK(key_press_event), &app);
    g_signal_connect(glarea, "realize", G_CALLBACK(gl_init), &app);
    g_signal_connect(glarea, "render", G_CALLBACK(gl_render), &app);

    g_signal_connect(window, "destroy", G_CALLBACK(close_window), &app);

    //fprintf(stderr, "calc surface\n");
    //calc_surface(&app);
    //fprintf(stderr, "calc surface done\n");

    gtk_widget_show_all (window);

    gtk_main();
    return 0;
}
