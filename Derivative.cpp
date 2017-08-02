#include <stdio.h>
#include "Derivative.h"


double Derivative2d_x(double **u,int i,int j,double h,char c)
{
    if(c == 'c')
        return (u[i+1][j]-u[i-1][j])/(2*h);
    else if(c == 'r')
        return (u[i+1][j]-u[i][j])/(h);
    else if(c == 'l')
        return (u[i][j]-u[i-1][j])/(h);
    else
        printf("Error Derivative x\n");
}
double Derivative2d_y(double **u,int i,int j,double h,char c)
{
    if(c == 'c')
        return (u[i][j+1]-u[i][j-1])/(2*h);
    else if(c == 'r')
        return (u[i][j+1]-u[i][j])/(h);
    else if(c == 'l')
        return (u[i][j]-u[i][j-1])/(h);
    else
        printf("Error Derivative y\n");
}
double Laplas2d_y(double **u,int i,int j,double h)
{
    return ((u[i][j+1]-2*u[i][j]+u[i][j-1])/(h*h));
}
double Laplas2d_x(double **u,int i,int j,double h)
{
    return ((u[i+1][j]-2*u[i][j]+u[i-1][j])/(h*h));
}
