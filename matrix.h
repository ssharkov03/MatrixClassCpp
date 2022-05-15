#ifndef MATRIX_H
#define MATRIX_H

#include <stdexcept>
#include <iostream>
#include <iomanip>

template <class T>
class Matrix
{
    public:
        // Define the various constructors
        Matrix();
        Matrix(int nRows, int nCols); // Initializes zero matrix of size nRows x nCols
        Matrix(int nRows, int nCols, const T *inputData);
        Matrix(const Matrix<T> &inputMatrix);

        // Destructor
        ~Matrix ();

        // Configuration methods
        void SetToIdentity();
        void SetToZero();

        // Element access methods
        T GetElement(int row, int col);
        T SetElement(int row, int col, T elementValue);
        int GetNumRows();
        int GetNumCols();

        // Overload == operator
        bool operator== (const Matrix<T> &right);

        // Overload +, -, * operators (friends)
        template<class U> friend Matrix<U> operator+ (const Matrix<U> &left, const Matrix<U> &right);
        template<class U> friend Matrix<U> operator+ (const U &left, const Matrix<U> &right);
        template<class U> friend Matrix<U> operator+ (const Matrix<U> &left, const U &right);

        template<class U> friend Matrix<U> operator- (const Matrix<U> &left, const Matrix<U> &right);
        template<class U> friend Matrix<U> operator- (const U &left, const Matrix<U> &right);
        template<class U> friend Matrix<U> operator- (const Matrix<U> &left, const U &right);

        template<class U> friend Matrix<U> operator* (const Matrix<U> &left, const Matrix<U> &right);
        template<class U> friend Matrix<U> operator* (const U &left, const Matrix<U> &right);
        template<class U> friend Matrix<U> operator* (const Matrix<U> &left, const U &right);

        // Overload = operator
        Matrix<T>& operator=(const Matrix& matrix);

        // Overload [] operator
        const T* operator[](int i) const;
        T& at(int i, int j);
        const T& at(int i, int j) const;

        // Matrix operations/transformations
        Matrix<T> Transpose();
        double Determinant();

        // Help functions
        void PrintMatrix();

    private:
        void BoundsCheck(int i, int j) const;
        bool IsSquare();
        int Sub2Ind(int row, int col);

    private:
        T *m_matrixData;
        int m_nRows, m_nCols, m_nElements;

};

#endif
/*---------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------*/




/*************************************************************************************
                            CONSTRUCTOR / DESTRUCTOR FUNCTIONS
*************************************************************************************/
// The default constructor
template <class T>
Matrix<T>::Matrix()
{
    m_nRows = 1;
    m_nCols = 1;
    m_nElements = 1;
    m_matrixData = new T[m_nElements];
    m_matrixData[0] = 0.0;
}


// Construct zero matrix
template <class T>
Matrix<T>::Matrix(int nRows, int nCols)
{
    m_nRows = nRows;
    m_nCols = nCols;
    m_nElements = nRows * nCols;
    m_matrixData = new T[m_nElements];
    for (int i = 0; i < m_nElements; ++i){
        m_matrixData[i] = 0.0;
    }
}


// Construct matrix from linear array
template <class T>
Matrix<T>::Matrix(int nRows, int nCols, const T *inputData)
{
    m_nRows = nRows;
    m_nCols = nCols;
    m_nElements = nRows * nCols;
    m_matrixData = new T[m_nElements];
    for (int i = 0; i < m_nElements; ++i){
        m_matrixData[i] = inputData[i];
    }
}

// The copy constructor
template <class T>
Matrix<T>::Matrix(const Matrix<T> &inputMatrix)
{
    m_nRows = inputMatrix.m_nRows;
    m_nCols = inputMatrix.m_nCols;
    m_nElements = inputMatrix.m_nElements;
    m_matrixData = new T[m_nElements];
    for (int i = 0; i < m_nElements; ++i){
        m_matrixData[i] = inputMatrix.m_matrixData[i];
    }
}

// Destructor
template <class T>
Matrix<T>::~Matrix()
{
    delete[] m_matrixData;
}

/*************************************************************************************
                            CONFIGURATION FUNCTIONS
*************************************************************************************/
// Eye matrix
template<class T>
void Matrix<T>::SetToIdentity(){

    if (!IsSquare()){
        throw std::invalid_argument("Cannot form an identity matrix that is not square.");
    }

    for(int row = 0; row < m_nRows; ++row){
        for(int col = 0; col < m_nCols; ++col){
            if (row == col){
                m_matrixData[Sub2Ind(row,col)] = 1.0;
            }
            else{
                m_matrixData[Sub2Ind(row,col)] = 0.0;
            }
        }
    }
}

// Zero matrix
template<class T>
void Matrix<T>::SetToZero(){

    for(int row = 0; row < m_nRows; ++row){
        for(int col = 0; col < m_nCols; ++col){
            m_matrixData[Sub2Ind(row,col)] = 0.0;
        }
    }
}

template<class T>
Matrix<T> Matrix<T>::Transpose() {
    int nRows = m_nRows;
    int nCols = m_nCols;
    Matrix<T> transposed_matrix(nCols, nRows);
    int m_LinearIdx, transposed_m_LinearIdx;

    // matrix[i][j] = transposed_matrix[j][i]
    for (int i = 0; i < nRows; ++i){
        for (int j = 0; j < nCols; ++j){
            m_LinearIdx = Sub2Ind(i,j);
            transposed_m_LinearIdx = (j * nRows) + i;  // Sub2Ind for transposed matrix
            transposed_matrix.m_matrixData[transposed_m_LinearIdx] = m_matrixData[m_LinearIdx];
        }
    }
    return transposed_matrix;
}

template<class T>
double Matrix<T>::Determinant() {
    if (!IsSquare()){
        throw std::invalid_argument("Cannot compute the determinant of matrix that is not square.");
    }
    else {

        // Copying original Matrix, so that it is not changed
        int sizeOfMatrix = m_nRows;
        Matrix<T> matrixForGaussMethod(sizeOfMatrix, sizeOfMatrix);

        for (int i = 0; i < sizeOfMatrix * sizeOfMatrix; ++i){
            matrixForGaussMethod.m_matrixData[i] = m_matrixData[i];
        }


        int fixedRow = 0, fixedColumnConsistsOfZeros, notZeroElemIdxInFixedColumn;
        double det = 1, firstElemOfFixedColumn, notZeroElemOfFixedColumn, transvectionAlpha;
        int linearIdx;

        for (int fixedColumn = 0; fixedColumn < sizeOfMatrix; ++fixedColumn) {

            fixedColumnConsistsOfZeros = 1;
            for (int row = fixedRow; row < sizeOfMatrix; ++row) {
                linearIdx = Sub2Ind(row, fixedColumn);
                if (matrixForGaussMethod.m_matrixData[linearIdx] != 0) {
                    fixedColumnConsistsOfZeros = 0;
                    break;
                }
            }
            if (fixedColumnConsistsOfZeros == 1) {
                continue;
            }
            linearIdx = Sub2Ind(fixedRow, fixedColumn);
            firstElemOfFixedColumn = matrixForGaussMethod.m_matrixData[linearIdx];
            notZeroElemIdxInFixedColumn = fixedColumn;
            if (firstElemOfFixedColumn == 0) {
                for (int row = fixedRow + 1; row < sizeOfMatrix; ++row) {
                    linearIdx = Sub2Ind(row, fixedColumn);
                    notZeroElemOfFixedColumn = matrixForGaussMethod.m_matrixData[linearIdx];
                    if (notZeroElemOfFixedColumn != 0) {
                        notZeroElemIdxInFixedColumn = row;
                    }
                }
                for (int column = 0; column < sizeOfMatrix; ++column) {
                    matrixForGaussMethod.m_matrixData[Sub2Ind(fixedRow, column)] += matrixForGaussMethod.m_matrixData[Sub2Ind(notZeroElemIdxInFixedColumn, column)];
                }
            }


            for (int row = fixedRow + 1; row < sizeOfMatrix; ++row) {
                transvectionAlpha = (-matrixForGaussMethod.m_matrixData[Sub2Ind(row, fixedColumn)]) / matrixForGaussMethod.m_matrixData[Sub2Ind(fixedRow, fixedColumn)];
                for (int column = fixedColumn; column < sizeOfMatrix; ++column) {
                    matrixForGaussMethod.m_matrixData[Sub2Ind(row, column)] = matrixForGaussMethod.m_matrixData[Sub2Ind(row, column)] + matrixForGaussMethod.m_matrixData[Sub2Ind(fixedRow, column)] * transvectionAlpha;
                }
            }

            ++fixedRow;
        }


        for (int diag_element = 0; diag_element < sizeOfMatrix; ++diag_element) {
            det *= matrixForGaussMethod.m_matrixData[Sub2Ind(diag_element, diag_element)];
        }

        return det;
    }
}


/*************************************************************************************
                            ELEMENT FUNCTIONS
*************************************************************************************/

template <class T>
T Matrix<T>::GetElement(int row, int col) {
    int linearIdx = Sub2Ind(row, col);
    if (linearIdx >= 0){
        return m_matrixData[linearIdx];
    }
    return 0.0;
}


template <class T>
T Matrix<T>::SetElement(int row, int col, T elementValue) {
    int linearIdx = Sub2Ind(row, col);
    if (linearIdx >= 0){
        m_matrixData[linearIdx] = elementValue;
        return true;
    }
    else {
        return false;
    }
}

template <class T>
int Matrix<T>::GetNumRows()
{
    return m_nRows;
}

template <class T>
int Matrix<T>::GetNumCols()
{
    return m_nCols;
}

/*************************************************************************************
                            OPERATOR FUNCTIONS
*************************************************************************************/

/*************************************************************************************
THE + OPERATOR
*************************************************************************************/

// matrix + matrix
template<class T>
Matrix<T> operator+ (const Matrix<T> &left, const Matrix<T> &right)
{
    int numRows = left.m_nRows;
    int numCols = left.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; ++i){
        tempResult[i] = left.m_matrixData[i] + right.m_matrixData[i];
    }
    Matrix<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// scalar + matrix
template<class T>
Matrix<T> operator+ (const T &left, const Matrix<T> &right)
{
    int numRows = right.m_nRows;
    int numCols = right.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; ++i){
        tempResult[i] = left + right.m_matrixData[i];
    }
    Matrix<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// matrix + scalar
template<class T>
Matrix<T> operator+ (const Matrix<T> &left, const T &right)
{
    int numRows = left.m_nRows;
    int numCols = left.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; ++i){
        tempResult[i] = left.m_matrixData[i] + right;
    }
    Matrix<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}


/*************************************************************************************
THE - OPERATOR
*************************************************************************************/

// matrix - matrix
template<class T>
Matrix<T> operator- (const Matrix<T> &left, const Matrix<T> &right)
{
    int numRows = left.m_nRows;
    int numCols = left.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; ++i){
        tempResult[i] = left.m_matrixData[i] - right.m_matrixData[i];
    }
    Matrix<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// scalar - matrix
template<class T>
Matrix<T> operator- (const T &left, const Matrix<T> &right)
{
    int numRows = right.m_nRows;
    int numCols = right.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; ++i){
        tempResult[i] = left - right.m_matrixData[i];
    }
    Matrix<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}


// matrix - scalar
template<class T>
Matrix<T> operator- (const Matrix<T> &left, const T &right)
{
    int numRows = left.m_nRows;
    int numCols = left.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; ++i){
        tempResult[i] = left.m_matrixData[i] - right;
    }
    Matrix<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

/*************************************************************************************
THE * OPERATOR
*************************************************************************************/

//scalar * matrix
template<class T>
Matrix<T> operator* (const T &left, const Matrix<T> &right)
{
    int numRows = right.m_nRows;
    int numCols = right.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; ++i){
        tempResult[i] = left * right.m_matrixData[i];
    }
    Matrix<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

//matrix * scalar
template<class T>
Matrix<T> operator* (const Matrix<T> &left, const T &right)
{
    int numRows = left.m_nRows;
    int numCols = left.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; ++i){
        tempResult[i] = left.m_matrixData[i] * right;
    }
    Matrix<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// matrix * matrix
template<class T>
Matrix<T> operator* (const Matrix<T> &left, const Matrix<T> &right)
{
    int l_numRows = left.m_nRows;
    int l_numCols = left.m_nCols;
    int r_numRows = right.m_nRows;
    int r_numCols = right.m_nCols;

    if (l_numCols == r_numRows){

        T *tempResult = new T[left.m_nRows * right.m_nCols];
        // fixating row in left and column in right
        for (int leftRow = 0; leftRow < l_numRows; ++leftRow){
            for (int rightCol = 0; rightCol < r_numCols; ++rightCol){
                T elementResult = 0.0;
                for (int leftCol = 0; leftCol < l_numCols; ++leftCol){  // loop through each elem in fixed vectors matrix
                    int leftLinearIdx = (leftRow * l_numCols) + leftCol;
                    int rightLinearIdx = (leftCol * r_numCols) + rightCol; // leftCol = rightRow
                    elementResult += (left.m_matrixData[leftLinearIdx] * right.m_matrixData[rightLinearIdx]);
                }
                int resultLinearIdx = (leftRow * r_numCols) + rightCol;
                tempResult[resultLinearIdx] = elementResult;
            }
        }
        Matrix <T> result(l_numRows, r_numCols, tempResult);
        delete[] tempResult;
        return result;
    }
    else{
        Matrix <T> result(1, 1);
        return result;
    }
}

/*************************************************************************************
THE == OPERATOR
*************************************************************************************/
template<class T>
bool Matrix<T>::operator== (const Matrix<T> &right) {

    // Check if matrices are the same size, if not return false.
    if ((this->m_nRows != right.m_nRows) && (this->m_nCols != right.m_nCols)){
        return false;
    }
    // Check if elements are equal
    for (int i = 0; i < this->m_nElements; ++i){
        if (this->m_matrixData[i] != right.m_matrixData[i]){
            return false;
        }
    }
    return true;
}

/*************************************************************************************
THE = OPERATOR
*************************************************************************************/


template<class T>
Matrix<T> & Matrix<T>::operator=(const Matrix<T> &matrix) {
    if (this == &matrix){
        return *this;
    }
    if (m_nRows != matrix.m_nRows || m_nCols != matrix.m_nCols){
        delete[] m_matrixData;
        m_nRows = matrix.m_nRows;
        m_nCols = matrix.m_nCols;
        m_nElements = m_nRows * m_nCols;
        m_matrixData = new T[m_nRows * m_nCols];
    }

    for (int i = 0; i < m_nElements; ++i){
        m_matrixData[i] = matrix.m_matrixData[i];
    }
    return *this;
}


/*************************************************************************************
THE [] OPERATOR
*************************************************************************************/

template<class T>
const T *Matrix<T>::operator[](int i) const {
    return &this->m_matrixData[i * m_nCols];
}

template<class T>
T &Matrix<T>::at(int i, int j) {
    BoundsCheck(i, j);
    return this->m_matrixData[i * m_nCols + j];
}

template<class T>
const T &Matrix<T>::at(int i, int j) const {
    BoundsCheck(i, j);
    return this->m_matrixData[i * m_nCols + j];
}



/*************************************************************************************
                            HELP FUNCTIONS
*************************************************************************************/

template <class T>
void Matrix<T>::PrintMatrix(){
    int nRows = this->GetNumRows();
    int nCols = this->GetNumCols();

    for (int row = 0; row < nRows; ++row){
        for (int col = 0; col < nCols; ++col){
            std::cout << std::fixed << std::setprecision(2) << this->GetElement(row, col) << " ";
        }
        std::cout << std::endl;
    }
}

/*************************************************************************************
                            PRIVATE FUNCTIONS
*************************************************************************************/

template<class T>
bool Matrix<T>::IsSquare() {
    if (m_nRows == m_nCols){
        return true;
    }
    else{
        return false;
    }
}

// For matrix[row][col] returns idx of element in linear array m_matrixData
template <class T>
int Matrix<T>::Sub2Ind(int row, int col)
{
    if ((row < m_nRows) && (row >= 0) && (col < m_nCols) && (col >= 0)){
        return (row * m_nCols) + col;
    }
    return -1;
}


template<class T>
void Matrix<T>::BoundsCheck(int i, int j) const {
    if (!(0 <= i && i < m_nRows && 0 <= j && j < m_nCols)) {
        throw std::out_of_range("Can not access out of bound element!");
    }
}















