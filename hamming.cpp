#include <iostream>
#include "hamming.h"



using namespace std;
int main(){

    int r=3, p=2;
    int n = (int) ((pow(p, r) - 1) / (p - 1));
    Hamming hamming_code(p,r);
    printf("\n--------Parity Check Matrix-------\n");
    hamming_code.ParityCheckMatrix().PrintMatrix();
    printf("\n--------Generator Matrix----------\n");
    hamming_code.GeneratorMatrix().PrintMatrix();
    printf("------------------------------------\n");

    int simpleData[4] = {1, 0, 1, 0};
    Matrix<int> u(1, 4, simpleData);
    Matrix<int> x = u * hamming_code.GeneratorMatrix();
    x.MatrixByMod(p);
    x.PrintMatrix();

    printf("\n-----Code word: y vector----------\n");
    Matrix<int> y = x;
    y.SetElement(0, 1, 1);
    y.PrintMatrix();

    printf("\n------Decoded y vector------------\n");
    hamming_code.decode(y).PrintMatrix();
    return 0;
}
