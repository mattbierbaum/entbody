#ifndef __PLOT_H__
#define __PLOT_H__

#include <GL/freeglut.h>

void plot_init(int size);
void plot_clean();

int *plot_render_particles(float *x, float *r, int *c, long N, float L, float *shade);
void plot_render_line(float L, float *p, float *v);
int plot_clear_screen();
int plot_exit_func();

void plot_init_opengl();
void plot_end_opengl();
void plot_set_draw_color(float r, float g, float b, float a);
 
#endif
