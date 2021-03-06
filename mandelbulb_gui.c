#include <gtk/gtk.h>
#include <math.h>
#include <limits.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

#include <epoxy/gl.h>

#include "mandelbulb.h"

struct App {
    GtkWidget* window;
    GtkWidget* glarea;
    GtkWidget* combo1;
    GtkWidget* combo2;
    int stride;
    const char* type;

    double tx, ty, tz;

    guint program;
    guint mvp_location;
    guint position_location;
    guint color_location;
    guint vao;

    float mvp[16];

    struct VertexInfo* vertex_data;
    long nvertex;
    struct Mandelbulb* mandelbulb;
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

static void calc_surface(struct App* app) {
    app->vertex_data = mandelbulb_surface_vertices(app->mandelbulb, &app->nvertex);
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
    fprintf(stderr, "nvertex %ld\n", app->nvertex);

    guint vao, buffer;
    guint color_index = app->color_location;
    guint position_index = app->position_location;
    struct VertexInfo* vertex_data = app->vertex_data;
    long nvertex = app->nvertex;

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

void cleanup(struct App* app) {
    if (app->vao) {
        glDeleteVertexArrays(1, &app->vao);
        app->vao = 0;
    }
    free(app->vertex_data); app->vertex_data = NULL;
    app->nvertex = 0;
}

void gl_cleanup(GtkWidget* widget, struct App* app)
{
    gtk_gl_area_make_current (GTK_GL_AREA(widget));
    cleanup(app);
    mandelbulb_free(app->mandelbulb); app->mandelbulb = NULL;
}

void combo_changed(GtkComboBox* widget, struct App* app) {
    if (
        strcmp(gtk_combo_box_get_active_id(GTK_COMBO_BOX(app->combo1)), app->type)
        || atoi(gtk_combo_box_get_active_id(GTK_COMBO_BOX(app->combo2))) != app->stride)
    {
        app->type = gtk_combo_box_get_active_id(GTK_COMBO_BOX(app->combo1));
        app->stride = atoi(gtk_combo_box_get_active_id(GTK_COMBO_BOX(app->combo2)));

        mandelbulb_free(app->mandelbulb); app->mandelbulb = NULL;
        if (!strcmp(app->type, "nine")) {
            app->mandelbulb = mandelbulb_nine_new(app->stride);
        } else if (!strcmp(app->type, "quintic")) {
            app->mandelbulb = mandelbulb_quintic_new(app->stride, 0, 0, 2, 0);
        } else if (!strcmp(app->type, "cubic")) {
            app->mandelbulb = mandelbulb_cubic_new(app->stride);
        } else if (!strcmp(app->type, "quadratic")) {
            app->mandelbulb = mandelbulb_quadratic_new(app->stride);
        }
    } else {
        return;
    }
    gtk_gl_area_make_current (GTK_GL_AREA(app->glarea));
    cleanup(app);
    init_buffers(app);
    gtk_widget_queue_draw(app->glarea);
}

static void close_window(GtkWidget* window, struct App* app)
{
    gtk_main_quit();
}

int main(int argc, char** argv) {
    struct App app;
    memset(&app, 0, sizeof(struct App));
    GtkWidget* window;
    GtkWidget* glarea;
    GtkWidget* box_main;
    GtkWidget* box_top;
    GtkWidget* combo1;
    GtkWidget* combo2;

    gtk_init(&argc, &argv);

    init_mvp (app.mvp);

    window = gtk_window_new (GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(window), "Mandelbulb");
    gtk_window_set_default_size(GTK_WINDOW(window), 800, 800);

    box_top = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
    box_main = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);

    glarea = gtk_gl_area_new();
    gtk_gl_area_set_has_depth_buffer(GTK_GL_AREA(glarea), TRUE);

    gtk_box_pack_end(GTK_BOX(box_main), glarea, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(box_main), box_top, FALSE, TRUE, 0);
    gtk_container_add(GTK_CONTAINER(window), box_main);

    combo1 = gtk_combo_box_text_new();
    combo2 = gtk_combo_box_text_new();
    gtk_box_pack_start(GTK_BOX(box_top), combo1, TRUE, TRUE, 0);
    gtk_box_pack_end(GTK_BOX(box_top), combo2, TRUE, TRUE, 0);

    gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(combo1), "quadratic", "quadratic formula");
    gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(combo1), "cubic", "cubic formula");
    gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(combo1), "quintic", "quintic formula");
    gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(combo1), "nine", "nine formula");
    gtk_combo_box_set_active_id(GTK_COMBO_BOX(combo1), "nine");

    gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(combo2), "128", "resolution 128x128x128");
    gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(combo2), "256", "resolution 256x256x256");
    gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(combo2), "512", "resolution 512x512x512");
    gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(combo2), "1024", "resolution 1024x1024x1024");
    gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(combo2), "1280", "resolution 1280x1280x1280");

    gtk_combo_box_set_active_id(GTK_COMBO_BOX(combo2), "256");
    app.stride = 256;
    app.mandelbulb = mandelbulb_nine_new(app.stride);
    app.type = "nine";

    app.combo1 = combo1;
    app.combo2 = combo2;
    app.window = window;
    app.glarea = glarea;

    g_signal_connect(combo1, "changed", G_CALLBACK(combo_changed), &app);
    g_signal_connect(combo2, "changed", G_CALLBACK(combo_changed), &app);
    g_signal_connect(window, "key-press-event", G_CALLBACK(key_press_event), &app);
    g_signal_connect(glarea, "realize", G_CALLBACK(gl_init), &app);
    g_signal_connect(glarea, "render", G_CALLBACK(gl_render), &app);

    g_signal_connect(glarea, "unrealize", G_CALLBACK(gl_cleanup), &app);
    g_signal_connect(window, "destroy", G_CALLBACK(close_window), &app);

    gtk_widget_show_all (window);

    gtk_main();
    return 0;
}
