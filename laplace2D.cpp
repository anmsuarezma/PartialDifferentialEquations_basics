#include <iostream>
#include <cstdlib>
#include <cmath>

const double lx = 2.;
const double ly = 2.;
const int nx = 30;
const int ny = 30; //Cantidad de puntos
const int nt = 400;
const double dx = lx/(nx-1);
const double dy = ly/(ny-1);

void evolucione(double u[nx][ny]);
void condicionesDeFrontera(double u[nx][ny]);
void imprimase(const double u[][ny]);
void condicionesIniciales(double u[nx][ny]);
void gnuplot(void);

int main(void)
{
    double p[nx][ny] = {};
    condicionesIniciales(p);
    condicionesDeFrontera(p);
    gnuplot();

    for (int n = 0; n <= nt; n++)
    {
        evolucione(p);
        //condicionesDeFrontera(p);
        imprimase(p);        
    }
    //imprimase(p);        
    return 0;
}

void condicionesDeFrontera(double p[nx][ny])
{
	
	/* Condiciones de Dirichlet */   
    //Izquierda
    for (int j = 0; j < ny; j++)
        p[0][j] = 0.0;       
        
    //Derecha
    for (int j = 0; j < ny; j++)
        p[nx-1][j] = 100.0;    
    

	/* Condiciones de Newmann*/
	for (int ii = 0; ii < nx; ++ii)
	{
		p[ii][1] = p[ii][0]; //y=1
		p[ii][0] = p[ii][ny-1]; //y=0
	}
        
}

void evolucione(double p[nx][ny])
{
    double newP[nx][ny] = {};
    double dx2 = pow(dx,2);
    double dy2 = pow(dy,2);

	 for (int i = 0; i < nx; ++i)
    	for (int j = 0; j < ny; j++)
        	newP[i][j] = p[i][j];           
           
        
	for (int i = 1; i < nx - 1 ; ++i)
    	for (int j = 1; j < ny - 1; j++)
            p[i][j] = (1. / (2. * (dx2 + dy2))) * (dx2 * (newP[i][j+1] + newP[i][j-1]) + dy2*(newP[i+1][j] + newP[i-1][j])); 
    
}
void imprimase(const double p[][ny])
{
    std::cout << "splot '-' w l lw 2 " << std::endl;
    double x, y;
    for (int ii = 0; ii < nx; ++ii)
    {
        x = ii * dx;
        for (int jj = 0; jj < ny; ++jj)
        {
            y = jj * dy;
            std::cout << x << "  " << y << "  " << p[ii][jj] << std::endl;
        }
        std::cout << std::endl;
    }
    std::cout << "e" << std::endl;
    //std::cout << "set zrange [1:2] " << std::endl;
}

void condicionesIniciales(double p[nx][ny])
{
    for (int i = 1; i < nx - 1; i++)
        for (int j = 1; j < ny - 1; ++j)
            p[i][j] =  - 50.0;              
        
    
}

void gnuplot(void)
{
    std::cout << "set terminal gif animate" << std::endl;
    std::cout << "set out 'laplace2D.gif'" << std::endl;
    std::cout << "set contour base" << std::endl;
    std::cout << "set pm3d" << std::endl;
 	std::cout << "set xlabel 'x'" << std::endl;
	std::cout << "set ylabel 'y'" << std::endl;
	std::cout << "set title 'Laplacian(Pressure)'" << std::endl;
}