#ifndef __FIELDS_H__
#define __FIELDS_H__

typedef struct {
    int size[3];
    int numcomps;
    float *data;
} field;

field *fields_calculate_scalar(int N, float L, float *x, float *s, int S);
field *fields_calculate_vector(int N, float L, float *x, float *v, int S);

float *fields_value_at(field *fld, float *x);
float *fields_dfdx_at(field *fld, float *x);
float *fields_dfdy_at(field *fld, float *x);
float *fields_grad_at(field *fld, float *x);
#endif
