#include <iostream>
#include <cstdlib>

const double dx = 0.1;
const double dy = 0.1;
const double dt = 0.01;
const double lx = 2.;
const double ly = 2.;
const int nx = lx / dx + 1;
const int ny = ly / dy + 1;
const int nt = 200;
const double c = 1.;

void diferenciasFinitas(double u[][ny]);
void condicionesDeFrontera(double u[][ny]);
void imprimase(const double u[][ny]);
void condicionesIniciales(double u[][ny]);
void gnuplot(void);

int main(void)
{
    double u[nx][ny] = {};
    condicionesIniciales(u);
    condicionesDeFrontera(u);
    gnuplot();

    for (int n = 0; n <= nt; n++)
    {
        diferenciasFinitas(u);
        condicionesDeFrontera(u);
        imprimase(u);
    }
    
    return 0;
}

void condicionesDeFrontera(double u[nx][ny])
{
    //Abajo
    for (int i=0; i<nx; i++)
        u[i][0] = 1.0;
    //Arriba
    for (int i=0; i<nx; i++)
        u[i][ny-1] = 1.0;
    //Izquierda
    for (int j=0; j<ny; j++)
        u[0][j] = 1.0;
     //Derecha
    for (int j=0; j<ny; j++)
        u[nx-1][j] = 1.0;        
}

void diferenciasFinitas(double u[nx][ny])
{
    double newU[nx][ny] = {};

    for (int i = 0; i < nx - 1; ++i)
        for (int j=0; j<ny -1;j++)
            newU[i][j] = u[i][j];
        


    for (int i = 0; i < nx - 1; ++i)
        for (int j=0; j<ny -1;j++)
            u[i][j] = -(c*dt/(dx*dx))*(newU[i][j] - newU[i-1][j]) -(c*dt/(dy*dy))*(newU[i][j+1] - newU[i][j-1])+ newU[i][j];
      
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

void condicionesIniciales(double u[nx][ny])
{
    for (int i=1; i<nx-1; i++)
    {
        for (int j=0; j<ny-1;++j)
        {
            u[i][j] = 1.0;
        }
    }
     for (int ii=nx/4; ii<nx/2; ii++)
    {
        for (int jj=ny/4; jj<ny/2;++jj)
        {
            u[ii][jj] = 2.0;
        }
    }
}

void gnuplot(void)
{
    std::cout << "set terminal gif animate" << std::endl;
    std::cout << "set out 'linearConvection2D.gif'" << std::endl;
    std::cout << "set contour base" << std::endl;
    std::cout << "set pm3d" << std::endl;
}