#include"pch.h"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <omp.h>
#include <fstream>
using namespace std;


double *u, **f; //переменные для параллельой реализации
double*u1, **f1;   //переменные для последовательной реализации
double eps = 0.0001;
double max;

//параллельная реализация закомментирована

void mat(int x)
{
	srand(time(0));
	f = new double *[x];
	       //f1 = new double *[x];
	for (int i = 0; i < x; ++i)
	{
		f[i] = new double[x + 1];
		    //f1[i] = new double[x + 1];
		f[i][x] = rand() % 15;
			//f[i][x]= f1[i][x] = rand() % 15;
	}
	int t = (x +4) / 2 - 1, dj = 0; 

	//t-переменная для сравнения с диагональным элементом
	//Меньше его на единицу

	for (int i = 0; i < x - 1; ++i)
	{
		dj = 0;
		for (int k = 0; k < i; ++k)
			dj += f[k][i];
		f[i][i] = (int)(x + 4) / 2;
			//f[i][i]= f1[i][i] = (int)(x + 4) / 2;
		for (int j = i + 1; j < x; ++j)
		{
			if (dj < t)
			{
				f[i][j] = f[j][i] = rand() % (t - dj);
					//f[i][j]= f1[i][j] = f1[j][i] = rand() % (t - dj);
				dj += f[i][j];
			}
			else
				f[i][j] = f[j][i]= 0;
				    //f[i][j] = f[j][i] = f1[i][j] = f1[j][i] = 0;
		}
	}
	f[x - 1][x - 1] = (int)(x + 4) / 2;
	    //f[x - 1][x - 1] = f1[x - 1][x - 1] = (int)(x + 4) / 2;
	u = new double[x];
	    //u1 = new double[x];
	for (int i = 0; i < x; ++i)
		u[i] = (rand() % 10)*0.07; //рандомно определяем начальные приближения
		 //u1[i] = u[i]= (rand()% 10)*0.07; //рандомно определяем начальные приближения

}
void OmpCalc(int l)
{
	max = 0;
	double t;
	// преобразование матрицы
	for (int i = 0; i < l; ++i)
	{
		t = -f[i][i];
		f[i][i] = 0;
		for (int j = 0; j < l + 1; ++j)
			f[i][j] = f[i][j] / t;
	}
	for (int i = 0; i < l; ++i)
		f[i][l] *= -1;

	// поиск решений
	do
	{
		max = 0;
#pragma omp parallel for num_threads(4) schedule(dynamic,10)
		for (int i = 0; i < l; i++)
		{
			double  u0 = 0;;
			for (int j = 0; j < l; j++)
				u0 += f[i][j] * u[j];
			u0 += f[i][l];
			double d = abs(u[i] - u0);
			u[i] = u0;
			if (d > max)
#pragma omp critical
				max = d;
		}
	} while (max > eps);
}
void Calc(int x)
{
	double t;
	for (int i = 0; i < x; ++i)
	{
		t = -f1[i][i];
		f1[i][i] = 0;
		for (int j = 0; j < x + 1; ++j)
			f1[i][j] = f1[i][j] / t;
	}
	for (int i = 0; i < x; ++i)
		f1[i][x] *= -1;
	max=0;
	do
	{
		max = 0;
		for (int i = 0; i < x; i++)
		{
			double u0 = 0;
			for (int j = 0; j < x; j++)
				u0 += f1[i][j] * u1[j];
			u0 += f1[i][x];
			double d = abs(u1[i] - u0);
			if (d > max)
				max = d;
			u1[i] = u0;
		}
	} while (max > eps);
}

double Checkparposl(int l)
{
	mat(l);
	OmpCalc(l);
	double temp = max;
	cout << "Checkparposl" << endl;
	Calc(l);
	return(abs(max - temp));
}
void Checkright(int l)
{
	mat(l);
	//Вывод сгенерированных матриц
	for (int i = 0; i < l; ++i)
	{
		for (int j = 0; j < l + 1; ++j)
			cout << f[i][j] << " ";
		cout << endl;
	}
	for (int i = 0; i < l; ++i)
	{
		for (int j = 0; j < l + 1; ++j)
			cout << f1[i][j] << " ";
		cout << endl;
	}

	// Решение
	OmpCalc(l);
	Calc(l);

	// Вывод решений
	for (int i = 0; i < l; ++i)
		cout << u1[i] << " ";
	cout << endl;

	for (int i = 0; i < l; ++i)
		cout << u[i] << " ";
	cout << endl;

}



int main()
{
	ofstream out1("C:\\tools\\55.txt");
	for (int l = 100; l <400; l += 250)
	{
		mat(l);


		double tt = omp_get_wtime();
		OmpCalc(l);
		tt = omp_get_wtime() - tt;

		out1.precision(10);
		out1 << fixed << tt << endl;

		for (int i = 0; i < 5; ++i)
			cout << u[i];
		delete[]u;     //Удаление выделенной памяти 
		for (int i = 0; i < l; ++i)
			delete[]f[i];
		delete[]f; 

		//Удаление выделенной для последовательной реализации памяти
		
		/*for (int i = 0; i < l; ++i)
			delete[]f1[i];
		delete[]f1;
		delete[]u1;*/
	}
}