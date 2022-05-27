#ifndef MATRIX_CPP_HAMMING_H
#define MATRIX_CPP_HAMMING_H

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "matrix.h"

class Hamming : public Matrix<int>{
    public:
        // Define constructor
        Hamming(int p, int r);  // p is alphabet {0, ... , p-1}, r is nRows in parity check matrix

        Matrix<int> ParityCheckMatrix();
        Matrix<int> GeneratorMatrix();
        Matrix<int> decode(const Matrix<int>& codeWordWithMistake);

    private:

        static Matrix<int> AddOneToVector(Matrix<int> vector, int vector_size, int mod);
        static bool isVectorCollinearToColumnsFromRange(const Matrix<int>& vector, const Matrix<int>& matrix, int idxBegin, int idxEnd, int mod);
        static bool isVectorCollinearToColumn(const Matrix<int>& vector, Matrix<int> matrix, int idx, int mod);
        static bool isCollinear(const Matrix<int>& vector1, const Matrix<int>& vector2, int vectorLen, int mod);
        static bool isCoefficientValid(int vector1_coord, int vector2_coord, double coefficient, bool coefficientWasMadeOfDiv_El2_On_El1, int mod);
        Matrix<int> findSyndrome(const Matrix<int>& parityCheckMatrix, Matrix<int> vector) const;
        static bool isSyndromeZeroColumn(Matrix<int> syndrome);
        static int findIdxWithMistake(Matrix<int> parityCheckMatrix, const Matrix<int>& syndrome);
        Matrix<int> findCodeWordWithoutMistakes(const Matrix<int>& parityCheckMatrix, const Matrix<int>& codeWordWithMistake);

    private:
        int P, R, N;
        Matrix<int> h_ParityCheckMatrix, h_GeneratorMatrix, leftSideParityCheckMatrix, rightSideGeneratorMatrix;
};

#endif
/*---------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/




/*************************************************************************************
                            CONSTRUCTOR FUNCTION
*************************************************************************************/

Hamming::Hamming(int p, int r) {
    P = p;
    R = r;
    N = (int) ((pow(P, R) - 1) / (P - 1));
    Matrix <int> zeroMatrix1(R, N);
    Matrix <int> zeroMatrix2(N - R, N);
    Matrix <int> zeroMatrix1left(R, N - R);
    Matrix <int> zeroMatrix2right(N - R, R);

    h_ParityCheckMatrix = zeroMatrix1;
    h_GeneratorMatrix = zeroMatrix2;
    leftSideParityCheckMatrix = zeroMatrix1left;
    rightSideGeneratorMatrix = zeroMatrix2right;

    int nRows = R;
    int nCols = N;

    Matrix <int> identityMatrix(nRows, nRows);
    identityMatrix.SetToIdentity();

    for (int j = nCols - nRows; j < nCols; ++j){
        for (int i = 0; i < nRows; ++i){
            h_ParityCheckMatrix.SetElement(i, j, identityMatrix[i][j - (nCols - nRows)]);
        }
    }
    Matrix <int> cur_vector(1, nRows);  // vector of R zeros
    cur_vector.SetElement(0, nRows-1, 1);

    for (int j = nCols - nRows - 1; j >= 0; --j){
        while (isVectorCollinearToColumnsFromRange(cur_vector, h_ParityCheckMatrix, j + 1, nCols, P)){
            cur_vector = AddOneToVector(cur_vector, R, P);
        }
        for (int i = 0; i < nRows; ++i){
            h_ParityCheckMatrix.SetElement(i, j, cur_vector[0][i]);
        }
    }

    int K = N - R;
    for (int i = 0; i < R; ++i){
        for (int j = 0; j < K; ++j){
            leftSideParityCheckMatrix.SetElement(i, j, h_ParityCheckMatrix[i][j]);
        }
    }
    rightSideGeneratorMatrix = leftSideParityCheckMatrix.Transpose();
    rightSideGeneratorMatrix = -1 * rightSideGeneratorMatrix;

    Matrix <int> identityMatrix2(K, K);
    identityMatrix2.SetToIdentity();


    for (int i = 0; i < K; ++i){
        for (int j = 0; j < K; ++j){
            h_GeneratorMatrix.SetElement(i, j, identityMatrix2[i][j]);
        }
        for (int j = K; j < N; ++j){
            h_GeneratorMatrix.SetElement(i, j, rightSideGeneratorMatrix[i][j - K]);
        }
    }
}

/*************************************************************************************
                            PUBLIC FUNCTIONS
*************************************************************************************/

Matrix<int> Hamming::ParityCheckMatrix() {
    return h_ParityCheckMatrix;
}


Matrix<int> Hamming::GeneratorMatrix() {
    return h_GeneratorMatrix;
}

Matrix<int> Hamming::decode(const Matrix<int>& codeWordWithMistake) {

    Matrix<int> codeWordWithoutMistake;
    codeWordWithoutMistake = findCodeWordWithoutMistakes(h_ParityCheckMatrix, codeWordWithMistake);

    // first n-r coords of correct coords are coords of decoded code word because generator matrix G has E(n-k*n-k) at beginning
    // and codeWord = decodedCodeWord * G

    Matrix<int> decodedCodeWord(1, N-R);
    for (int i = 0; i < N - R; ++i){
        decodedCodeWord.SetElement(0, i, codeWordWithoutMistake[0][i]);
    }
    return decodedCodeWord;
}


/*************************************************************************************
                            PRIVATE FUNCTIONS
*************************************************************************************/

Matrix<int> Hamming::AddOneToVector(Matrix<int> vector, int vector_size, int mod) {
    for (int i = vector_size - 1; i >= 0; --i){
        if (vector[0][i] + 1 != mod){
            vector.SetElement(0, i, vector[0][i] + 1);
            return vector;
        }
        else{
            if (i == 0){
                std::cout << "Could not increase vector!" << std::endl;
                exit(-1);
            }
            vector.SetElement(0, i, 0);
        }
    }
    return vector;
}


bool Hamming::isVectorCollinearToColumnsFromRange(const Matrix<int>& vector, const Matrix<int>& matrix, int idxBegin, int idxEnd, int mod) {
    for (int idx = idxBegin; idx < idxEnd; ++idx){
        if (isVectorCollinearToColumn(vector, matrix, idx, mod)){
            return true;
        }
    }
    return false;
}


bool Hamming::isVectorCollinearToColumn(const Matrix<int>& vector, Matrix<int> matrix, int idx, int mod) {

    int vectorLen = matrix.GetNumRows();
    Matrix<int> columnRemake(1, vectorLen);
    for (int i = 0; i < vectorLen; ++i){
        columnRemake.SetElement(0, i, matrix[i][idx]);
    }
    if (isCollinear(vector, columnRemake, vectorLen, mod)){
        return true;
    }
    return false;
}


bool Hamming::isCollinear(const Matrix<int>& vector1, const Matrix<int>& vector2, int vectorLen, int mod) {
    double coefficient = NULL;
    bool coefficientWasMadeOfDiv_El2_On_El1 = NULL;

    for (int i = 0; i < vectorLen; ++i){

        if ((vector1[0][i] + vector2[0][i] != 0) && (vector1[0][i] * vector2[0][i] == 0)){
            return false;
        }
        else if (vector1[0][i] == vector2[0][i] && vector2[0][i] == 0){
            continue;
        }
        else if (coefficient == NULL){
            if (vector2[0][i] > vector1[0][i]){
                coefficient = (double) vector2[0][i] / vector1[0][i];
                coefficientWasMadeOfDiv_El2_On_El1 = true;
            }
            else{
                coefficient = (double) vector1[0][i] / vector2[0][i];
                coefficientWasMadeOfDiv_El2_On_El1 = false;
            }
        }
        else if (!isCoefficientValid(vector1[0][i], vector2[0][i], coefficient, coefficientWasMadeOfDiv_El2_On_El1, mod)){
            return false;
        }
    }
    return true;
}


bool Hamming::isCoefficientValid(int vector1_coord, int vector2_coord, double coefficient, bool coefficientWasMadeOfDiv_El2_On_El1, int mod) {
    if (coefficientWasMadeOfDiv_El2_On_El1){
        if (((int)(vector1_coord * coefficient) % mod == vector2_coord) && (vector1_coord * coefficient - (int)vector1_coord * coefficient == 0.0))
            return true;
        return false;
    }
    else{
        if (((int)(vector2_coord * coefficient) % mod == vector1_coord) && (vector2_coord * coefficient - (int)vector2_coord * coefficient == 0.0))
            return true;
        return false;
    }
}



Matrix<int> Hamming::findSyndrome(const Matrix<int>& parityCheckMatrix, Matrix<int> vector) const {
    Matrix<int> transposed_vector, syndrome;
    transposed_vector = vector.Transpose();
    syndrome = (parityCheckMatrix * transposed_vector);
    syndrome.MatrixByMod(P);
    return syndrome;
}

bool Hamming::isSyndromeZeroColumn(Matrix<int> syndrome) {

    // nRows = r, nCols = 1
    int nRows = syndrome.GetNumRows();
    for (int i = 0; i < nRows; ++i){
        if (syndrome[i][0] != 0){
            return false;
        }
    }
    return true;
}



int Hamming::findIdxWithMistake(Matrix<int> parityCheckMatrix, const Matrix<int>& syndrome) {

    // syndrome is a column in parity check matrix, here its idx is found -> mistake in idx coord of code word y(vector)

    int nCols = parityCheckMatrix.GetNumCols();
    int nRows = parityCheckMatrix.GetNumRows();

    for (int columnIdx = 0; columnIdx < nCols; ++columnIdx){
        bool isColumnEqToSyndrome = true;

        for (int rowIdx = 0; rowIdx < nRows; ++rowIdx){
            if (syndrome[rowIdx][0] != parityCheckMatrix[rowIdx][columnIdx]){
                isColumnEqToSyndrome = false;
                break;
            }
        }
        if (isColumnEqToSyndrome) {
            return columnIdx;
        }
    }
    return -1;
}

Matrix<int> Hamming::findCodeWordWithoutMistakes(const Matrix<int>& parityCheckMatrix, const Matrix<int>& codeWordWithMistake) {

    Matrix<int> syndrome;
    syndrome = findSyndrome(parityCheckMatrix, codeWordWithMistake);
    if (isSyndromeZeroColumn(syndrome)){
        return codeWordWithMistake;
    }

    int mistakeIdx = findIdxWithMistake(parityCheckMatrix, syndrome), mistakeIdxCoordValue = -1;
    Matrix<int> codeWordWithoutMistake;
    codeWordWithoutMistake = codeWordWithMistake;
    codeWordWithoutMistake.SetElement(0, mistakeIdx, mistakeIdxCoordValue);  // mistake idx search start (from 0)


    do{
        ++mistakeIdxCoordValue;
        codeWordWithoutMistake.SetElement(0, mistakeIdx, mistakeIdxCoordValue);
        syndrome = findSyndrome(parityCheckMatrix, codeWordWithoutMistake);
    }
    while (!isSyndromeZeroColumn(syndrome));

    // syndrome = 0 => codeWordWithoutMistake is correct now;
    return codeWordWithoutMistake;
}
