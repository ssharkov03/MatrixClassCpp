#include <iostream>
#include "hamming.h"

using namespace std;
int main(){

    int r=4, p=3;
    Hamming hamming_code(p,r);
    printf("\n--------Parity Check Matrix-------\n");
    hamming_code.ParityCheckMatrix().PrintMatrix();
    printf("\n--------Generator Matrix----------\n");
    hamming_code.GeneratorMatrix().PrintMatrix();

    return 0;
}