/**
 * Mandelbrot and Julia sets
 * Copyright (c) 2021, Alexey Ozeritskiy
 */

#include <gtk/gtk.h>
#include <math.h>
#include <stdlib.h>
#include "colormap_vga1.h"
#include "big_float.h"

using T = BigFloat<2,uint64_t>;
//using T = double;

struct App {
    GdkPixbuf* pixbuf;
    GtkWidget* window;
    GtkWidget* entry;
    GtkWidget* drawing_area;
    T julia_x;
    T julia_y;
    double t;
    T zoom;
    T zoom_initial;
    T zoom_x;
    T zoom_y;
    int zoom_xx;
    int zoom_yy;
    int julia;
    int max_iteration;
    int timer_id;
    int dirty;
};

template<typename T>
double get_iteration_mandelbrot(T x0, T y0, struct App* app) {
    T x = 0;
    T y = 0;
    T xn;
    int i;
    double fract_i;
    for (i = 1; i<app->max_iteration && x*x+y*y < 4; i=i+1) {
        if constexpr(std::is_same_v<T,double>) {
            xn = x*x - y*y + x0;
            y = T(2.0)*x*y + y0;
            x = xn;
        } else {
            xn = x*x - y*y + x0;
            y = (x*y).Mul2() + y0;
            x = xn;
        }
    }

    fract_i = i;
    if (i < app->max_iteration) {
        if constexpr(std::is_same_v<T, double>) {
            double log_zn = log((x*x + y*y)) / 2;
            double nu = log(log_zn / log(2)) / log(2);
            fract_i = fract_i + 1 - nu;
        } else {
            double log_zn = log((x*x + y*y).ToDouble()) / 2;
            double nu = log(log_zn / log(2)) / log(2);
            fract_i = fract_i + 1 - nu;
        }
    }

    return fract_i;
}

template<typename T>
double get_iteration_julia(T x0, T y0, struct App* app) {
    int i;
    T cx = app->julia_x;
    T cy = app->julia_y;
    T x = x0;
    T y = y0;
    T xn;
    T yn;
    double fract_i;
    for (i = 1; i<app->max_iteration && x*x+y*y < 1<<16; i=i+1) {
        xn = x*x - y*y + cx;
        yn = T(2.0)*x*y + cy;
        x = xn; y = yn;
    }

    fract_i = i;
    if (i < app->max_iteration) {
        if constexpr(std::is_same_v<T, double>) {
            double log_zn = log((x*x + y*y)) / 2;
            double nu = log(log_zn / log(2)) / log(2);
            fract_i = fract_i + 1 - nu;
        } else {
            double log_zn = log((x*x + y*y).ToDouble()) / 2;
            double nu = log(log_zn / log(2)) / log(2);
            fract_i = fract_i + 1 - nu;
        }
    }

    return fract_i;
}

static void from_screen_(
    T* x0, T* y0, T screen_x, T screen_y, int width, int height, T zoom, T zoom_x, T zoom_y) {
    T x = screen_x;
    T y = screen_y;
    int w = width < height ? width : height; /* keep aspect ratio */
    T _1_w = 1.0/w;
    T d = abs(width - height);
    /* scaling */
    x = x*_1_w;
    y = y*_1_w;
    d = d*_1_w;
    x = x+T(-0.5);
    y = y+T(-0.5);
    x = x*T(4);
    y = y*T(4);
    d = d*T(4);
    /* centering */
    if (width > height) {
        x = x+T(-0.5)*d;
    } else {
        y = y+T(-0.5)*d;
    }
    x = (x + (-zoom_x)) * zoom + zoom_x;
    y = (y + (-zoom_y)) * zoom + zoom_y;
    // TODO: fixed zoom point
    *x0 = x; *y0 = y;
}

static void from_screen(
    T* x0, T* y0, T screen_x, T screen_y, int width, int height, T zoom, T zoom_x, T zoom_y, int zoom_xx, int zoom_yy)
{
    T x1, y1;
    T x2, y2;
    from_screen_(&x1, &y1, screen_x, screen_y, width, height, zoom, 0, 0);
    from_screen_(&x2, &y2, zoom_xx, zoom_yy, width, height, zoom, 0, 0);
    *x0 = x1 - x2 + zoom_x;
    *y0 = y1 - y2 + zoom_y;
}

static void hsv2rgb(guchar* rgb, int h, double s, double v) {
    double r[3] = {0};
    int h1 = h/60; h1 %= 6;
    double vmin = (1-s)*v;
    double a = (v-vmin)*(h%60)/60.0;
    double vinc = vmin+a;
    double vdec = v-a;

    switch (h1) {
        case 0:
            r[0] = v; r[1] = vinc; r[2] = vmin;
            break;
        case 1:
            r[0] = vdec; r[1] = v; r[2] = vmin;
            break;
        case 2:
            r[0] = vmin; r[1] = v; r[2] = vinc;
            break;
        case 3:
            r[0] = vmin; r[1] = vdec; r[2] = v;
            break;
        case 4:
            r[0] = vinc; r[1] = vmin; r[2] = v;
            break;
        case 5:
            r[0] = v; r[1] = vmin; r[2] = vdec;
            break;
        default:
            break;
    }

    for (int i = 0; i < 3; i=i+1) {
        rgb[i] = (guchar)(r[i]*255);
    }

   // printf("%d %f %f -> %d %d %d,  %f %f %f\n", h, s, v, rgb[0], rgb[1], rgb[2], r[0], r[1], r[2]);
}

static void linear_interpolate(guchar* p, double it, int max_iteration) {
    int off = 150;
    double fractpart, intpart;
    it /= max_iteration;
    it *= 360;
    fractpart = modf(it, &intpart);
    int j = round(intpart);
    int h1 = 20*j;
    int h2 = 20*(j+1);
    int h = h1 + fractpart * (h2 - h1);
    hsv2rgb(p, off+h1, 1, 1);
    h = off+h;
    if (h > 360) {
        h %= 360;
    }
    hsv2rgb(p, h, 1, 1);
}

template<typename T>
void print_zoom(T zoom) {
    if constexpr(std::is_same_v<T, double>) {
        std::cerr << "zoom=" << zoom << std::endl;
    } else {
        std::cerr << "zoom=" << zoom.ToDouble() << std::endl;
    }
}

static void create_pixbuf(GtkWidget *widget, struct App* app)
{
    if (app->pixbuf) {
        g_object_unref(app->pixbuf);
    }

    int width = gtk_widget_get_allocated_width (widget);
    int height = gtk_widget_get_allocated_height (widget);

    if (app->zoom_xx < 0) {
        app->zoom_xx = width / 2;
        app->zoom_yy = height / 2;
    }

    GdkPixbuf *	pixbuf = gdk_pixbuf_new(GDK_COLORSPACE_RGB, FALSE, 8, width, height);

    int n_channels = gdk_pixbuf_get_n_channels (pixbuf);
    int rowstride = gdk_pixbuf_get_rowstride (pixbuf);
    guchar* pixels = gdk_pixbuf_get_pixels (pixbuf);

    double (*get_iteration)(T x0, T y0, struct App* app);
    if (app->julia) {
        get_iteration = get_iteration_julia;
    } else {
        get_iteration = get_iteration_mandelbrot;
    }

    print_zoom(app->zoom);
#pragma omp parallel for
    for (int y = 0; y < height; y=y+1) {
        for (int x = 0; x < width; x=x+1) {
            T x0, y0;
            from_screen(&x0, &y0, x, y, width, height, app->zoom, app->zoom_x, app->zoom_y, app->zoom_xx, app->zoom_yy);

            guchar* p = pixels + y * rowstride + x * n_channels;
            double it = get_iteration(x0, y0, app);

            if (it < app->max_iteration) {
                // linear_interpolate(p, it, app->max_iteration);
                p[0] = 0;
                p[1] = 255;
                p[2] = 0;
            } else {
                p[0] = 255; // r
                p[1] = 0;
                p[2] = 0;
            }
        }
    }
    app->pixbuf = pixbuf;
}

static gboolean configure_event_cb(GtkWidget *widget, GdkEventConfigure *event, struct App* app)
{
    create_pixbuf(widget, app);
    return TRUE;
}

template<typename T>
void print_points(char* buf, T x, T y) {
    if constexpr(std::is_same_v<T, double>) {
        snprintf(buf, sizeof(buf), "%.4le + %.4le i", x, y);
    } else {
        snprintf(buf, sizeof(buf), "%.4le + %.4le i", x.ToDouble(), y.ToDouble());
    }
}

static gboolean motion_notify_event_cb(GtkWidget *widget, GdkEventMotion *event, struct App* app)
{
    int width = gtk_widget_get_allocated_width (widget);
    int height = gtk_widget_get_allocated_height (widget);
    T x, y;
    from_screen(&x, &y, event->x, event->y, width, height, app->zoom, app->zoom_x, app->zoom_y, app->zoom_xx, app->zoom_yy);
    char buf[1024];

    // fixed point
    app->zoom_x = x;
    app->zoom_y = y;
    app->zoom_xx = event->x;
    app->zoom_yy = event->y;

    print_points(buf, x, y);
    gtk_entry_set_text(GTK_ENTRY(app->entry), buf);

    return TRUE;
}

static gboolean button_press_event_cb(GtkWidget *widget, GdkEventButton *event, struct App* app)
{
    double x, y;
    char title[1024];
    const char* text = gtk_entry_get_text(GTK_ENTRY(app->entry));
    sscanf(text, "%le + %le i", &x, &y);
    sprintf(title, "Julia: %s", text);
    app->julia_x = x;
    app->julia_y = y;
    app->julia = 1;
    create_pixbuf(GTK_WIDGET (app->drawing_area), app);
    gtk_widget_queue_draw (GTK_WIDGET (app->drawing_area));
    gtk_window_set_title(GTK_WINDOW(app->window), title);
    return TRUE;
}

static gboolean reset_button_released_event_cb(GtkWidget *btn, GdkEvent  *event, struct App* app)
{
    if (app->timer_id > 0)
    {
        g_source_remove(app->timer_id);
        app->timer_id = 0;
    }
    app->t = 0;
    app->julia = 0;
    app->zoom = 1.0;
    app->zoom_x = 0.0;
    app->zoom_y = 0.0;
    app->zoom_xx = -1;
    create_pixbuf(GTK_WIDGET (app->drawing_area), app);
    gtk_widget_queue_draw (GTK_WIDGET (app->drawing_area));
    gtk_window_set_title(GTK_WINDOW(app->window), "Mandelbrot");
    return TRUE;
}

static gboolean redraw_timeout(struct App *app)
{
    app->t += 0.001;
    app->julia_x = 2*cos(7*app->t);
    app->julia_y = 2*sin(11*app->t);
    create_pixbuf(GTK_WIDGET (app->drawing_area), app);
    char buf[1024];

    print_points(buf, app->julia_x, app->julia_y);
    gtk_entry_set_text(GTK_ENTRY(app->entry), buf);

    char title[2048];
    sprintf(title, "Julia: %s", buf);
    gtk_window_set_title(GTK_WINDOW(app->window), title);

    if (!app->dirty) {
        gtk_widget_queue_draw (GTK_WIDGET (app->drawing_area));
    }
    app->dirty = 1;
    return TRUE;
}

static gboolean run_button_released_event_cb(GtkWidget *btn, GdkEvent* event, struct App* app)
{
    app->julia = 1;
    app->t = 0;
    app->julia_x = 0;
    app->julia_y = 0;
    if (app->timer_id > 0)
    {
        g_source_remove(app->timer_id);
        app->timer_id = 0;
    }
    app->timer_id = g_timeout_add(100, (GSourceFunc)redraw_timeout, app);
    return TRUE;
}

static gboolean draw_cb(GtkWidget *widget, cairo_t *cr, struct App*  app)
{
    gdk_cairo_set_source_pixbuf(cr, app->pixbuf, 0, 0);
    cairo_paint(cr);
    app->dirty = 0;

    return TRUE;
}

static void
zoom_begin_cb (GtkGesture       *gesture,
               GdkEventSequence *sequence,
               struct App         *app)
{
    app->zoom_initial = app->zoom;
}

static void
zoom_end_cb (GtkGesture       *gesture,
             GdkEventSequence *sequence,
             struct App         *app)
{
}

static void
zoom_scale_changed_cb (GtkGestureZoom *z,
                       gdouble         scale,
                       struct App       *app)
{
    app->zoom = app->zoom_initial * T(1.0 / scale);
    create_pixbuf(GTK_WIDGET (app->drawing_area), app);
    if (!app->dirty) {
        gtk_widget_queue_draw (GTK_WIDGET (app->drawing_area));
    }
    app->dirty = 1;
}

static gboolean mouse_scroll (GtkWidget *widget,
               GdkEvent  *event,
               struct App       *app)
{
    GdkEventScroll * scroll = (GdkEventScroll*) event;

    if (scroll->direction) {
        app->zoom = app->zoom * 1.5;
    } else {
        app->zoom = app->zoom * 0.5;
    }
    create_pixbuf(GTK_WIDGET (app->drawing_area), app);
    gtk_widget_queue_draw (GTK_WIDGET (app->drawing_area));

    return TRUE;
}

static void close_window(GtkWidget* widget, struct App* app)
{
    if (app->timer_id > 0)
    {
        g_source_remove(app->timer_id);
        app->timer_id = 0;
    }

    if (app->pixbuf)
    {
        g_object_unref(app->pixbuf);
    }

    gtk_main_quit();
}

int main(int argc, char** argv) {
    struct App app;
    memset(&app, 0, sizeof(struct App));
    GtkWidget* window;
    GtkWidget* drawing_area;
    GtkWidget* box; /*main box*/
    GtkWidget* box_right; /*right box with button*/
    GtkWidget* box_buttons;
    GtkWidget* entry;
    GtkWidget* button;
    GtkWidget* button2;
    GtkGesture * zoom;

    gtk_init(&argc, &argv);

    window = gtk_window_new (GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(window), "Mandelbrot");
    gtk_window_set_default_size(GTK_WINDOW(window), 800, 600);
    drawing_area = gtk_drawing_area_new();
    zoom = gtk_gesture_zoom_new(drawing_area);
    box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
    box_right = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
    box_buttons = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
    entry = gtk_entry_new();
    button = gtk_button_new_with_label("reset");
    button2 = gtk_button_new_with_label("run");

    gtk_container_add(GTK_CONTAINER(window), box);
    gtk_box_pack_start(GTK_BOX(box), drawing_area, TRUE, TRUE, 0);
    gtk_box_pack_end(GTK_BOX(box), box_right, FALSE, TRUE, 0);

    gtk_box_pack_start(GTK_BOX(box_right), entry, FALSE, TRUE, 0);
    gtk_box_pack_end(GTK_BOX(box_right), box_buttons, FALSE, TRUE, 0);

    gtk_box_pack_start(GTK_BOX(box_buttons), button, FALSE, TRUE, 0);
    gtk_box_pack_end(GTK_BOX(box_buttons), button2, FALSE, TRUE, 0);

    gtk_box_set_child_packing(GTK_BOX(box_right), box_buttons, FALSE, TRUE, 0, GTK_PACK_START);

    gtk_widget_set_events(drawing_area,
                          gtk_widget_get_events(drawing_area)
                          | GDK_BUTTON_PRESS_MASK | GDK_SCROLL_MASK | GDK_POINTER_MOTION_MASK);

    app.window = window;
    app.entry = entry;
    app.drawing_area = drawing_area;
    app.max_iteration = 1000;
    app.zoom = 1.0;
    app.zoom_xx = app.zoom_yy = -1;
    app.zoom_initial = 1.0;

    /* readonly entry */
    GValue val = G_VALUE_INIT;
    g_value_init(&val, G_TYPE_BOOLEAN);
    g_value_set_boolean(&val, FALSE);
    g_object_set_property(G_OBJECT(entry), "editable", &val);
    g_object_set_property(G_OBJECT(entry), "can_focus", &val);

    g_signal_connect(drawing_area, "draw", G_CALLBACK(draw_cb), &app);
    g_signal_connect(drawing_area, "configure-event", G_CALLBACK(configure_event_cb), &app);
    g_signal_connect(drawing_area, "motion-notify-event", G_CALLBACK(motion_notify_event_cb), &app);
    g_signal_connect(drawing_area, "button-press-event", G_CALLBACK(button_press_event_cb), &app);
    g_signal_connect(drawing_area, "scroll-event", G_CALLBACK(mouse_scroll), &app); // zoom
    g_signal_connect(button, "button-release-event", G_CALLBACK(reset_button_released_event_cb), &app);
    g_signal_connect(button2, "button-release-event", G_CALLBACK(run_button_released_event_cb), &app);

    g_signal_connect(window, "destroy", G_CALLBACK(close_window), &app);

    g_signal_connect(zoom, "begin", G_CALLBACK(zoom_begin_cb), &app);
    g_signal_connect(zoom, "scale-changed", G_CALLBACK(zoom_scale_changed_cb), &app);
    g_signal_connect(zoom, "end", G_CALLBACK(zoom_end_cb), &app);

    gtk_widget_show_all (window);

    gtk_main();
    return 0;
}
