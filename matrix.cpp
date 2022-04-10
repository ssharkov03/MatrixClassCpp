#include <iostream>
#include "matrix.h"

using namespace std;
int main(){

    cout << "Code to test Matrix Class" << endl;

    double simpleData[12] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    Matrix <double> testMatrix(3, 4, simpleData);
    cout << endl << "**************" << endl;
    cout << "3x4 matrix test (testMatrix)" << endl;
    testMatrix.PrintMatrix();


    cout << endl << "**************" << endl;
    cout << "Test element retrieval" << endl;
    cout << "Element (0,0) = " << testMatrix.GetElement(0,0) << endl;
    cout << "Element (1,0) = " << testMatrix.GetElement(1,0) << endl;
    cout << "Element (2,0) = " << testMatrix.GetElement(2,0) << endl;
    cout << "Element (0,1) = " << testMatrix.GetElement(0,1) << endl;
    cout << "Element (1,1) = " << testMatrix.GetElement(1,1) << endl;
    cout << "Element (2,1) = " << testMatrix.GetElement(2,1) << endl;
    cout << "Element (5,5) = " << testMatrix.GetElement(5,5) << endl;


    cout << endl << "**************" << endl;
    cout << "Test matrix multiplication" << endl;
    double simpleData2[12] = {1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3};
    Matrix <double> testMatrix2(4, 3, simpleData2);
    cout << "4x3 matrix test (testMatrix2)" << endl;
    testMatrix2.PrintMatrix();
    cout << "Multiplication (testMatrix * testMatrix2) result: " << endl;
    Matrix <double> multiTest1 = testMatrix * testMatrix2;
    multiTest1.PrintMatrix();

    cout << endl << "**************" << endl;
    cout << "testMatrix2 * testMatrix: " << endl;
    Matrix <double> multiTest2 = testMatrix2 * testMatrix;
    multiTest2.PrintMatrix();

    cout << endl << "**************" << endl;
    cout << "Test multiplication of column vector by matrix" << endl;
    double columnData[3] = {1.5, 2.5, 3.5};
    double squareData[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
    Matrix<double> testColumn(3, 1, columnData);
    Matrix<double> squareMatrix(3, 3, squareData);
    cout << "Column vector =" << endl;
    testColumn.PrintMatrix();
    cout << "Square matrix =" << endl;
    squareMatrix.PrintMatrix();
    cout << "Column vector * Square matrix =" << endl;
    operator*(testColumn, squareMatrix).PrintMatrix();
    cout << endl;
    (testColumn * squareMatrix).PrintMatrix();
    cout << "Square matrix * Column vector =" << endl;
    (squareMatrix * testColumn).PrintMatrix();
    cout << "Square matrix + 1 =" << endl;
    (squareMatrix + 1.0).PrintMatrix();

    cout << endl << "**************" << endl;
    Matrix<int> testForEye(5,5);
    cout << "Test for setting to identity matrix" << endl;
    cout << "Test for eye = " << endl;
    testForEye.PrintMatrix();
    testForEye.SetToIdentity();
    cout << "Identity matrix = " << endl;
    testForEye.PrintMatrix();

    cout << endl << "**************" << endl;
    cout << "Test for = operation" << endl;
    int lhsArray[9] = {1,2,3,4,5,6,7,8,9};
    int rhsArray[9] = {1,1,1,1,1,1,1,1,1};
    Matrix<int> lhs(3,3,lhsArray);
    Matrix<int> rhs(3,3,rhsArray);
    Matrix<int> result;
    result = lhs + rhs;
    cout << "Left matrix = " << endl;
    lhs.PrintMatrix();
    cout << "Right matrix = " << endl;
    rhs.PrintMatrix();
    cout << "Left + Right matrix = " << endl;
    result.PrintMatrix();





    cout << endl << "**************" << endl;
    cout << "Test for transposing matrix " << endl;
    int dataForTranspose[10] = {1,2,3,4,5,6,7,8,9,10};
    Matrix<int> matrixForTranspose(2, 5, dataForTranspose);
    Matrix<int> transposedMatrix;
    transposedMatrix = matrixForTranspose.Transpose();
    cout << "Matrix to transpose = " << endl;
    matrixForTranspose.PrintMatrix();
    cout << "Transposed matrix = " << endl;
    transposedMatrix.PrintMatrix();
    cout << "Test if matrix to transpose doesnt change:" << endl;
    matrixForTranspose.PrintMatrix();
    cout << endl;
    cout << "Second test" << endl;
    int dataForTranspose2[49];
    for (int i = 1; i <= 49; ++i){
        dataForTranspose2[i - 1] = i;
    }
    Matrix<int> matrixForTranspose2(7, 7, dataForTranspose2);
    Matrix<int> transposedMatrix2;
    transposedMatrix2 = matrixForTranspose2.Transpose();
    cout << "Matrix to transpose 2 = " << endl;
    matrixForTranspose2.PrintMatrix();
    cout << "Transposed matrix 2 = " << endl;
    transposedMatrix2.PrintMatrix();

    cout << endl << "**************" << endl;
    cout << "Test for determinant of matrix " << endl;
    double detData[4] = {1,2,3,4};
    Matrix <double> matrixForDeterminant(2,2,detData);
    cout << "Matrix to calc det from = " << endl;
    matrixForDeterminant.PrintMatrix();
    double det = matrixForDeterminant.Determinant();
    cout << "Matrix determinant = " << det << endl;

    cout << endl << "**************" << endl;
    cout << "Test for [] operator overloading " << endl;
    cout << "Matrix = " << endl;
    matrixForDeterminant.PrintMatrix();
    cout << "Matrix[1][1] =  " << matrixForDeterminant[1][1] << endl;


    return 0;
}