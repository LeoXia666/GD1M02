#include <iostream>
#include <vector>
#include <fstream>
#include <stdexcept>

using namespace std;

class Matrix {
public:
    vector<vector<double>> data;
    int rows, cols;

    Matrix(int r, int c) : rows(r), cols(c) {
        data.resize(r, vector<double>(c, 0));
    }

    void display() const {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                cout << data[i][j] << " ";
            }
            cout << endl;
        }
    }
};

// Helper function to compute the determinant of a 3x3 matrix
double determinant3x3(const Matrix& A) {
    return A.data[0][0] * (A.data[1][1] * A.data[2][2] - A.data[1][2] * A.data[2][1]) -
        A.data[0][1] * (A.data[1][0] * A.data[2][2] - A.data[1][2] * A.data[2][0]) +
        A.data[0][2] * (A.data[1][0] * A.data[2][1] - A.data[1][1] * A.data[2][0]);
}

// Function to compute the determinant of a 4x4 matrix using Laplace expansion
double determinant(const Matrix& A) {
    if (A.rows != 4 || A.cols != 4) {
        throw runtime_error("Determinant implemented only for 4x4 matrices.");
    }

    double det = 0;
    for (int p = 0; p < 4; p++) {
        Matrix subMat(3, 3);
        for (int i = 1; i < 4; i++) {
            int subCol = 0;
            for (int j = 0; j < 4; j++) {
                if (j == p) continue;
                subMat.data[i - 1][subCol] = A.data[i][j];
                subCol++;
            }
        }
        det += (p % 2 == 0 ? 1 : -1) * A.data[0][p] * determinant3x3(subMat);
    }
    return det;
}

// Function to compute the cofactor of a matrix element
Matrix cofactorMatrix(const Matrix& A) {
    Matrix cof(A.rows, A.cols);
    for (int row = 0; row < 4; row++) {
        for (int col = 0; col < 4; col++) {
            Matrix subMat(3, 3);
            int subRow = 0;
            for (int i = 0; i < 4; i++) {
                if (i == row) continue;
                int subCol = 0;
                for (int j = 0; j < 4; j++) {
                    if (j == col) continue;
                    subMat.data[subRow][subCol] = A.data[i][j];
                    subCol++;
                }
                subRow++;
            }
            cof.data[row][col] = ((row + col) % 2 == 0) ? determinant3x3(subMat) : -determinant3x3(subMat);
        }
    }
    return cof;
}

// Function to transpose a matrix
Matrix transpose(const Matrix& A) {
    Matrix result(A.cols, A.rows);
    for (int i = 0; i < A.rows; i++) {
        for (int j = 0; j < A.cols; j++) {
            result.data[j][i] = A.data[i][j];
        }
    }
    return result;
}

// Function to compute the inverse of a 4x4 matrix
Matrix inverse(const Matrix& A) {
    if (A.rows != 4 || A.cols != 4) {
        throw runtime_error("Inverse implemented only for 4x4 matrices.");
    }

    double det = determinant(A);
    if (det == 0) {
        throw runtime_error("Matrix is singular and cannot be inverted.");
    }

    // Find the cofactor matrix
    Matrix cofMat = cofactorMatrix(A);

    // Transpose the cofactor matrix to get the adjugate matrix
    Matrix adjugate = transpose(cofMat);

    // Divide adjugate by determinant to get the inverse matrix
    Matrix inv(A.rows, A.cols);
    for (int i = 0; i < A.rows; i++) {
        for (int j = 0; j < A.cols; j++) {
            inv.data[i][j] = adjugate.data[i][j] / det;
        }
    }

    return inv;
}

// Scalar multiplication
Matrix scalarMultiply(const Matrix& A, double scalar) {
    Matrix result(A.rows, A.cols);
    for (int i = 0; i < A.rows; i++) {
        for (int j = 0; j < A.cols; j++) {
            result.data[i][j] = A.data[i][j] * scalar;
        }
    }
    return result;
}

// Matrix addition
Matrix add(const Matrix& A, const Matrix& B) {
    Matrix result(A.rows, A.cols);
    for (int i = 0; i < A.rows; i++) {
        for (int j = 0; j < A.cols; j++) {
            result.data[i][j] = A.data[i][j] + B.data[i][j];
        }
    }
    return result;
}

// Matrix subtraction
Matrix subtract(const Matrix& A, const Matrix& B) {
    Matrix result(A.rows, A.cols);
    for (int i = 0; i < A.rows; i++) {
        for (int j = 0; j < A.cols; j++) {
            result.data[i][j] = A.data[i][j] - B.data[i][j];
        }
    }
    return result;
}

// Matrix multiplication
Matrix multiply(const Matrix& A, const Matrix& B) {
    Matrix result(A.rows, B.cols);
    for (int i = 0; i < A.rows; i++) {
        for (int j = 0; j < B.cols; j++) {
            result.data[i][j] = 0;
            for (int k = 0; k < A.cols; k++) {
                result.data[i][j] += A.data[i][k] * B.data[k][j];
            }
        }
    }
    return result;
}

// Generate identity matrix
Matrix identityMatrix(int size) {
    Matrix result(size, size);
    for (int i = 0; i < size; i++) {
        result.data[i][i] = 1;
    }
    return result;
}

// Read matrices and scalar from file
void readMatrixFromFile(const string& filename, Matrix& A, Matrix& B, double& scalar) {
    ifstream infile(filename);
    if (!infile.is_open()) {
        throw runtime_error("Cannot open file.");
    }

    // Set hardcoded dimensions for matrices A and B (assuming 4x4 in your case)
    A.rows = 4;
    A.cols = 4;
    B.rows = 4;
    B.cols = 4;

    // Resize the matrices based on the known dimensions
    A.data.resize(A.rows, vector<double>(A.cols));
    B.data.resize(B.rows, vector<double>(B.cols));

    // Read the values for matrix A
    for (int i = 0; i < A.rows; i++) {
        for (int j = 0; j < A.cols; j++) {
            infile >> A.data[i][j];
        }
    }

    // Read the values for matrix B
    for (int i = 0; i < B.rows; i++) {
        for (int j = 0; j < B.cols; j++) {
            infile >> B.data[i][j];
        }
    }

    // Read the scalar value
    infile >> scalar;

    infile.close();
}



int main() {
    Matrix A(4, 4);  // Adjusted to 4x4 matrices
    Matrix B(4, 4);
    double scalar;
    string filename;

    cout << "Enter the filename to read matrices and scalar value: ";
    cin >> filename;
    // Reading matrices and scalar from the file
    readMatrixFromFile(filename, A, B, scalar);

    cout << "Matrix A:" << endl;
    A.display();
    cout << "Matrix B:" << endl;
    B.display();
    cout << "Scalar value: " << scalar << endl;

    // 1. |A|: Determinant of matrix A
    cout << "1) Determinant of A: " << determinant(A) << endl;

    // 2. Aᵀ: Transpose of matrix A
    cout << "2) Transpose of A:" << endl;
    Matrix result = transpose(A);
    result.display();

    // 3. A⁻¹: Inverse of matrix A
    try {
        cout << "3) Inverse of A:" << endl;
        result = inverse(A);
        result.display();
    }
    catch (const exception& e) {
        cerr << e.what() << endl;
    }

    // 4. Multiply A by scalar
    cout << "4) A multiplied by scalar " << scalar << ":" << endl;
    result = scalarMultiply(A, scalar);
    result.display();

    // 5. A + B: Matrix addition
    cout << "5) A + B:" << endl;
    result = add(A, B);
    result.display();

    // 6. A - B: Matrix subtraction
    cout << "6) A - B:" << endl;
    result = subtract(A, B);
    result.display();

    // 7. A * B: Matrix multiplication
    cout << "7) A * B:" << endl;
    result = multiply(A, B);
    result.display();

    // 8. B * A: Matrix multiplication
    cout << "8) B * A:" << endl;
    result = multiply(B, A);
    result.display();

    // 9. Generate an identity matrix
    cout << "9) Identity matrix (4x4):" << endl;
    result = identityMatrix(4);
    result.display();

    return 0;
}
