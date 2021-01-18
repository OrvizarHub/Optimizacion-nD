#include <iostream>
#include <vector>

#include "matrix.h"

using namespace std;

Matrix::Matrix()
{
    
}

Matrix::Matrix(int rows, int cols)
{
    vector<vector<double>> aux(rows, vector<double>(cols, 0));
    data = aux;
}

Matrix::Matrix(vector<vector<double>> aux)
{
    data = aux;
}

int Matrix::Rows()
{
    return data.size();
}

int Matrix::Cols()
{
    return data[0].size();
}

double Matrix::GetValue(int i, int j)
{
    return data[i][j];
}

void Matrix::SetValue(int i, int j, double aux)
{
    data[i][j] = aux;
}

void Matrix::AddValue(vector<double> aux)
{
    data.push_back(aux);
}

Matrix Matrix::GetRow(int i)
{
    int cols = data[i].size();
    Matrix aux(1, cols);

    for (int j = 0; j < cols; j++)
    {
        aux.SetValue(0,j,data[i][j]);
    }
    
    return aux;
}

void Matrix::SetRow(int i, Matrix x)
{
    for (int j = 0; j < x.Cols(); j++)
    {
        data[i][j] = x.GetValue(0,j);
    }   
}

Matrix Matrix::GetCol(int i)
{
    int rows = data.size();
    Matrix aux(rows, 1);

    for (int j = 0; j < rows; j++)
    {
        aux.SetValue(j,0,data[i][j]);
    }

    return aux;
}

void Matrix::SetCol(int i, Matrix x)
{
    for (int j = 0; j < x.Rows(); j++)
    {
        data[i][j] = x.GetValue(j,0);
    }
}