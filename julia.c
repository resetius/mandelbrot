/**
 * Mandelbrot and Julia sets 
 * Copyright (c) 2021, Alexey Ozeritskiy
 */

#include <gtk/gtk.h>
#include <math.h>

struct JuliaParams {
    double x;
    double y;
};

struct JuliaSurface {
    cairo_surface_t *surface;
    struct JuliaParams params;
};

struct App {
    struct JuliaSurface* surfaces;
    int nsurfaces;
};

static int cols = 4;
static int rows = 8;

/* параметры из книги С. Деменок "Фрактал: между мифом и ремеслом" */
static struct JuliaParams params[] = {
    {-0.39054, -0.58679}, {-0.07, -0.7}, {-0.747, 0.107}, {0.318623, 0.044899},
    {0.268545, -0.003483}, {-1.028382, -0.264756}, {-1.258842, 0.065330}, {0.27856, -0.003483},
    {0.138341, 0.649857}, {-0.392488, -0.587966}, {-0.192175, 0.656734}, {0.238498, 0.519198},
    {0.345, -0.365}, {0.350, -0.360}, {0.352, -0.362}, {0.32, 0.043}, 
    {0, 0}, {-0.25, 0}, {-0.5, 0}, {-0.7, 0},
    {-0.75, 0}, {-0.8, 0}, {-0.9, 0}, {-1, 0},
    {-0.24, 0.74}, {-0.12, 0.74}, {-0.012, 0.74}, {-0.0012, 0.74},
    {-0.11031, -1}, {-0.11, -0.63}, {-0.11, -0.635}, {0.108294, -0.670487}
};

static void clear_surface(cairo_surface_t *surface)
{
    cairo_t *cr;
    cr = cairo_create(surface);

    cairo_set_source_rgb(cr, 1, 1, 1);
    cairo_paint(cr);

    cairo_destroy(cr);
}

static gboolean configure_event_cb(GtkWidget *widget, GdkEventConfigure *event, struct JuliaSurface* j)
{
    if (j->surface)
        cairo_surface_destroy(j->surface);

    j->surface = gdk_window_create_similar_surface(gtk_widget_get_window(widget), CAIRO_CONTENT_COLOR,
                                                gtk_widget_get_allocated_width(widget),
                                                gtk_widget_get_allocated_height(widget));

    clear_surface(j->surface);

    return TRUE;
}

int get_iteration_julia(double x0, double y0, struct JuliaSurface* j) {
    int i;
    double cx = j->params.x;
    double cy = j->params.y;
    double x = x0;
    double y = y0;
    double xn;
    double yn;
    for (i = 1; i<32; i=i+1) {
        xn = x*x - y*y + cx;
        yn = 2*x*y + cy;
        if (sqrt(xn*xn+yn*yn) > 2) {
            return i;
        }
        x = xn; y = yn;
    }
    return 0;
}

static void from_screen(double* x0, double* y0, double screen_x, double screen_y, int width, int height) {
    double x = screen_x;
    double y = screen_y;
    int w = width < height ? width : height; /* keep aspect ratio */
    double d = abs(width - height);
    /* scaling */
    x /= w; y /= w; d /= w;
    x -= 0.5; y -= 0.5;
    x *= 4; y *= 4; d *= 4;
    y = -y; /* mirror */
    /* centering */
    if (width > height) {
        x -= 0.5*d;
    } else {
        y -= 0.5*d;
    }
    *x0 = x; *y0 = y;
}

static gboolean draw_cb(GtkWidget *widget, cairo_t *cr, struct JuliaSurface* j)
{
    double x0, y0;

    int width = gtk_widget_get_allocated_width (widget);
    int height = gtk_widget_get_allocated_height (widget);

    GdkPixbuf *	pixbuf = gdk_pixbuf_new(GDK_COLORSPACE_RGB, FALSE, 8, width, height);

    int n_channels = gdk_pixbuf_get_n_channels (pixbuf);
    int rowstride = gdk_pixbuf_get_rowstride (pixbuf);
    int x, y;
    guchar* pixels = gdk_pixbuf_get_pixels (pixbuf);
    guchar* p;

    int it;

    for (y = 0; y < height; y=y+1) {
        for (x = 0; x < width; x=x+1) {
            from_screen(&x0, &y0, x, y, width, height);
 
            p = pixels + y * rowstride + x * n_channels;
            it = get_iteration_julia(x0, y0, j);

            if (it == 0) {
                p[0] = 255; // r
                p[1] = 0;
                p[2] = 0;
            } else {
                p[0] = it*8;
                p[1] = it*8;
                p[2] = (32-it)*8;
            }
        }
    }

    gdk_cairo_set_source_pixbuf(cr, pixbuf, 0, 0);

    g_object_unref(pixbuf);

    cairo_paint(cr);

    return FALSE;
}

static void close_window(GtkWidget *widget, struct App* app)
{
    int i;
    for (i = 0; i < app->nsurfaces; i=i+1) {
        if (app->surfaces[i].surface) {
            cairo_surface_destroy(app->surfaces[i].surface);
        }
    }

    gtk_main_quit();
}

int main(int argc, char** argv) {
    int row;
    int col;
    struct App app;
    GtkWidget* window;
    GtkWidget* grid;
    gtk_init(&argc, &argv);

    window = gtk_window_new (GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(window), "Julia");
    gtk_window_set_default_size(GTK_WINDOW(window), 600, 400);
    grid = gtk_grid_new();
    gtk_grid_set_row_homogeneous(GTK_GRID(grid), TRUE);
    gtk_grid_set_column_homogeneous(GTK_GRID(grid), TRUE);
    gtk_container_add(GTK_CONTAINER(window), grid);

    app.surfaces = malloc(rows*cols*sizeof(struct JuliaSurface));
    for (row = 0; row < rows; row=row+1) {
        for (col = 0; col < cols; col=col+1) {
            char buf[256];
            struct JuliaSurface* j = &app.surfaces[row*cols + col];
            j->surface = NULL;
            j->params = params[row*cols + col];
            if (j->params.y >= 0) {
                snprintf(buf, 256, "%.6lf+%.6lf i", j->params.x, j->params.y);
            } else {
                snprintf(buf, 256, "%.6lf%.6lf i", j->params.x, j->params.y);
            }
            GtkWidget* box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
            GtkWidget* drawing_area = gtk_drawing_area_new();
            GtkWidget* label = gtk_label_new(buf);

            gtk_box_pack_start(GTK_BOX(box), label, FALSE, TRUE, 0);
            gtk_box_pack_end(GTK_BOX(box), drawing_area, TRUE, TRUE, 0);

            gtk_grid_attach(GTK_GRID(grid), box, col, row, 1, 1);

            g_signal_connect(drawing_area, "draw", G_CALLBACK(draw_cb), j);
            g_signal_connect(drawing_area, "configure-event", G_CALLBACK(configure_event_cb), j);
        }
    }

    app.nsurfaces = rows*cols;
    g_signal_connect(window, "destroy", G_CALLBACK(close_window), &app);

    gtk_widget_show_all (window);

    gtk_main();
    return 0;
}
