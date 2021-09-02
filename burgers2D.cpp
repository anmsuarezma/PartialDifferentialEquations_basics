#include <iostream>
#include <cstdlib>
#include <cmath>

const double dx = 0.1;
const double dy = 0.1;
const double dt = 0.01;
const double lx = 2.;
const double ly = 2.;
const int nx = lx / dx + 1;
const int ny = ly / dy + 1;
const int nt = 120;
const double nu = 0.01;

void diferenciasFinitas(double u[nx][ny], double v[nx][ny]);
void condicionesDeFrontera(double u[nx][ny], double v[nx][ny]);
void imprimase(const double u[][ny]);
void condicionesIniciales(double u[nx][ny], double v[nx][ny]);
void gnuplot(void);

int main(void)
{
    double u[nx][ny] = {};
    double v[nx][ny] = {};
    condicionesIniciales(u,v);
    condicionesDeFrontera(u,v);
    gnuplot();

    for (int n = 0; n <= nt; n++)
    {
        diferenciasFinitas(u,v);
        condicionesDeFrontera(u, v);
        imprimase(u);
    }
    return 0;
}

void condicionesDeFrontera(double u[nx][ny], double v[nx][ny])
{
    //Abajo
    for (int i=0; i<nx; i++)
        {
            u[i][0] = 1.0;        v[i][0] = 1.0;
        }
    //Arriba
    for (int i=0; i<nx; i++)
        {
            u[i][ny-1] = 1.0;    v[i][ny-1] = 1.0;
        }
    //Izquierda
    for (int j=0; j<ny; j++)
        {
            u[0][j] = 1.0;       v[0][j] = 1.0;
        }
    //Derecha
    for (int j=0; j<ny; j++)
        {
            u[nx-1][j] = 1.0;    v[nx-1][j] = 1.0;
        }
}

void diferenciasFinitas(double u[nx][ny], double v[nx][ny])
{
    double newU[nx][ny] = {};
    double newV[nx][ny] = {};

    for (int i=0; i<nx; ++i)
    {
        for (int j=0; j<ny; ++j)
        {
            newU[i][j] = u[i][j];         
            newV[i][j] = v[i][j];
        }
    }

    for (int i=0; i<nx; ++i)
    {
        for (int j=0; j<ny; ++j)
        {
            u[i][j] = newU[i][j] + (nu*dt/pow(dx,2))*(newU[i+1][j]-2*newU[i][j]+newU[i-1][j]) + (nu*dt/pow(dy,2))*(newU[i][j+1]-2*newU[i][j]+newU[i][j-1]) - (dt/dx)*newU[i][j]*(newU[i][j]-newU[i-1][j]) - (dt/dy)*newV[i][j]*(newU[i][j]-newU[i][j-1]); 
            v[i][j] = newV[i][j] + (nu*dt/pow(dx,2))*(newV[i+1][j]-2*newV[i][j]+newV[i-1][j]) + (nu*dt/pow(dy,2))*(newV[i][j+1]-2*newV[i][j]+newV[i][j-1]) - (dt/dx)*newU[i][j]*(newV[i][j]-newV[i-1][j]) - (dt/dy)*newV[i][j]*(newV[i][j]-newV[i][j-1]); 
        }
    }
}
void imprimase(const double u[][ny])
{
        
    std::cout << "splot '-' w l lw 2 " << std::endl;

    double x,y;
    for (int ii=0; ii<nx; ++ii)
    {
        x = ii * dx;
        for (int jj=0; jj<ny; ++jj)
        {
            y = jj*dy;
            std::cout << x << "  " << y << "  " << u[ii][jj] << std::endl;
        }
        std::cout << std::endl;
    }

    std::cout << "e" << std::endl;
    std::cout << "set zrange [1:2] " << std::endl;
}

void condicionesIniciales(double u[nx][ny], double v[nx][ny])
{
    for (int i=1; i<nx-1; i++)
    {
        for (int j=0; j<ny-1;++j)
        {
            u[i][j] = 1.0;              v[i][j] = 1.0;
        }
    }
     for (int ii=nx/4; ii<nx/2+1; ii++)
    {
        for (int jj=ny/4; jj<ny/2+1;++jj)
        {
            u[ii][jj] = 2.0;            u[ii][jj] = 2.0;
        }
    }
}

void gnuplot(void)
{
    std::cout << "set terminal gif animate" << std::endl;
    std::cout << "set out 'burgers2D.gif'" << std::endl;
    std::cout << "set contour base" << std::endl;
    std::cout << "set pm3d" << std::endl;
}