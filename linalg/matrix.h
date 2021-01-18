#pragma once 

#include <vector>

using namespace std;

class Matrix
{
    private:
        vector<vector<double>> data;
    
    public:
        Matrix();
        Matrix(int rows, int cols);
        Matrix(vector<vector<double>> aux);

        int Rows();
        int Cols();

        double GetValue(int i, int j);
        void SetValue(int i, int j, double aux);
        void AddValue(vector<double>);
        Matrix GetRow(int i);
        void SetRow(int i, Matrix x);
        Matrix GetCol(int i);
        void SetCol(int i, Matrix x);
};