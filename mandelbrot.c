/**
 * Mandelbrot and Julia sets 
 * Copyright (c) 2021, Alexey Ozeritskiy
 */

#include <gtk/gtk.h>
#include <math.h>
#include <stdlib.h>
#include "colormap_vga1.h"

struct App {
    GdkPixbuf* pixbuf;
    GtkWidget* window;
    GtkWidget* entry;
    GtkWidget* drawing_area;
    double julia_x;
    double julia_y;
    int julia;
};

int get_iteration_mandelbrot(double x0, double y0, struct App* app) {
    double x = 0;
    double y = 0;
    double xn, yn;
    int i;
    for (i = 1; i<32; i=i+1) {
        xn = x*x - y*y + x0;
        yn = 2*x*y + y0;
        if (xn*xn+yn*yn > 4) {
            return i;
        }
        x = xn; y = yn;
    }
    return 0;
}

int get_iteration_julia(double x0, double y0, struct App* app) {
    int i;
    double cx = app->julia_x;
    double cy = app->julia_y;
    double x = x0;
    double y = y0;
    double xn;
    double yn;
    for (i = 1; i<200; i=i+1) {
        xn = x*x - y*y + cx;
        yn = 2*x*y + cy;
        if (xn*xn+yn*yn > 4) {
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
    /* centering */
    if (width > height) {
        x -= 0.5*d;
    } else {
        y -= 0.5*d;
    }
    *x0 = x; *y0 = y;
}

static void create_pixbuf(GtkWidget *widget, struct App* app)
{
    if (app->pixbuf) {
        g_object_unref(app->pixbuf);
    }

    int width = gtk_widget_get_allocated_width (widget);
    int height = gtk_widget_get_allocated_height (widget);

    GdkPixbuf *	pixbuf = gdk_pixbuf_new(GDK_COLORSPACE_RGB, FALSE, 8, width, height);

    int n_channels = gdk_pixbuf_get_n_channels (pixbuf);
    int rowstride = gdk_pixbuf_get_rowstride (pixbuf);
    int x, y;
    guchar* pixels = gdk_pixbuf_get_pixels (pixbuf);
    guchar* p;

    int it;
    double x0, y0;

    int (*get_iteration)(double x0, double y0, struct App* app);
    if (app->julia) {
        get_iteration = get_iteration_julia;
    } else {
        get_iteration = get_iteration_mandelbrot;
    }
    for (y = 0; y < height; y=y+1) {
        for (x = 0; x < width; x=x+1) {
            from_screen(&x0, &y0, x, y, width, height);
 
            p = pixels + y * rowstride + x * n_channels;
            it = get_iteration(x0, y0, app);
            if (it == 0) {
                p[0] = 255; // r
                p[1] = 0;
                p[2] = 0;
            } else {
                int off = 50;
                p[0] = colormap_vga1[it+off][0];
                p[1] = colormap_vga1[it+off][1];
                p[2] = colormap_vga1[it+off][2];
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

static gboolean motion_notify_event_cb(GtkWidget *widget, GdkEventMotion *event, struct App* app)
{
    int width = gtk_widget_get_allocated_width (widget);
    int height = gtk_widget_get_allocated_height (widget);
    double x, y;
    from_screen(&x, &y, event->x, event->y, width, height);
    char buf[1024];
    snprintf(buf, sizeof(buf), "%.4le + %.4le i", x, y);
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

static gboolean reset_button_released_event_cb(GtkButton *btn, struct App* app)
{
    app->julia = 0;
    gtk_widget_queue_draw (GTK_WIDGET (app->drawing_area));
    gtk_window_set_title(GTK_WINDOW(app->window), "Mandelbrot");
    return TRUE;
}

static gboolean draw_cb(GtkWidget *widget, cairo_t *cr, struct App*  app)
{
    gdk_cairo_set_source_pixbuf(cr, app->pixbuf, 0, 0);
    cairo_paint(cr);

    return TRUE;
}

static void close_window(GtkWidget* widget, struct App* app)
{   
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
    GtkWidget* entry;
    GtkWidget* button;
    gtk_init(&argc, &argv);

    window = gtk_window_new (GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(window), "Mandelbrot");
    gtk_window_set_default_size(GTK_WINDOW(window), 600, 400);
    drawing_area = gtk_drawing_area_new();
    box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
    box_right = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
    entry = gtk_entry_new();
    button = gtk_button_new_with_label("reset");

    gtk_container_add(GTK_CONTAINER(window), box);
    gtk_box_pack_start(GTK_BOX(box), drawing_area, TRUE, TRUE, 0);
    gtk_box_pack_end(GTK_BOX(box), box_right, FALSE, TRUE, 0);

    gtk_box_pack_start(GTK_BOX(box_right), entry, FALSE, TRUE, 0);

    gtk_box_pack_end(GTK_BOX(box_right), button, FALSE, TRUE, 0);
    gtk_box_set_child_packing(GTK_BOX(box_right), button, FALSE, TRUE, 0, GTK_PACK_START);

    gtk_widget_set_events(drawing_area,
                          gtk_widget_get_events(drawing_area) 
                          | GDK_BUTTON_PRESS_MASK | GDK_SCROLL_MASK | GDK_POINTER_MOTION_MASK);

    app.window = window;
    app.entry = entry;
    app.drawing_area = drawing_area;
    
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
    g_signal_connect(button, "released", G_CALLBACK(reset_button_released_event_cb), &app);

    g_signal_connect(window, "destroy", G_CALLBACK(close_window), &app);

    gtk_widget_show_all (window);

    gtk_main();
    return 0;
}
