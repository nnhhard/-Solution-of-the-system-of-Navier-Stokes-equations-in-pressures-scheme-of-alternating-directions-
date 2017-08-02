#ifndef DERIVATIVE_H_INCLUDED
#define DERIVATIVE_H_INCLUDED

double Derivative2d_x(double **u,int i,int j,double h,char c);
double Derivative2d_y(double **u,int i,int j,double h,char c);

double Laplas2d_x(double **u,int i,int j,double h);
double Laplas2d_y(double **u,int i,int j,double h);

#endif // DERIVATIVE_H_INCLUDED
