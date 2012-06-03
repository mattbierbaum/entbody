#ifndef __FIELDS_H__
#define __FIELDS_H__

typedef struct {
    int size[2];
    int numcomps;
    float *data;
} field;

field *fields_calculate_density(int N, float L, float *x);
field *fields_calculate_velocity(int N, float L, float *x, float *v);

float *fields_value_at(field *fld, float *x);
float *fields_dfdx_at(field *fld, float *x);
float *fields_dfdy_at(field *fld, float *x);
float *fields_grad_at(field *fld, float *x);
#endif
