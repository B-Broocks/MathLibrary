#include <iostream>
#include <vector>
#include <stdexcept>

using std::cout;
using std::endl;
using std::vector;

class Matrix {
private:
    int rows;
    int cols;
    vector<vector<double>> data; //matrix data
public:
    Matrix(int r, int c) : rows(r), cols(c) {
        //Initialize matrix with zero first
        data.resize(rows, vector<double>(cols, 0.0));
    }
    // plus
    Matrix operator+(const Matrix& a){
        if(rows != a.rows || cols != a.cols){
            throw std::invalid_argument("Matrix dimension mismatch!");
        }
        Matrix result(rows, cols);
        for(int i = 0; i < a.rows; ++i){
            for(int j = 0; j < a.cols; ++j){
                result.data[i][j] = data[i][j] + a.data[i][j];
            }
        }
        return result;
    }
    // Minus
    Matrix operator-(const Matrix& a){
        if(rows != a.rows || cols != a.cols){
            throw std::invalid_argument("Matrix dimension mismatch!");
        }
        Matrix result(rows, cols);
        for(int i = 0; i < a.rows; ++i){
            for(int j = 0; j < a.cols; ++j){
                result.data[i][j] = data[i][j] - a.data[i][j];
            }
        }
        return result;
    }
    // Mal
    Matrix operator*(double a){
        Matrix result(rows, cols);
        for(int i = 0; i < rows; ++i){
            for(int j = 0; j < cols; ++j){
                result.data[i][j] = data[i][j] * a;
            }
        }
        return result;
    }

    Matrix operator*(const Matrix& other) const{
        if(cols != other.rows){
            throw std::invalid_argument("Matrix dimensions incompatible for multiplication!");
        }
        Matrix result(rows, other.cols);

        for(int i = 0; i < rows; ++i){
            for(int j = 0; j < other.cols; ++j){

                double sum = 0;
                for(int k = 0; k < cols; ++k){

                    sum += data[i][k] * other.data[k][j];
                }
                result.data[i][j] = sum;
            }
        }
        return result;
    }



    Matrix identity_matrix() const{
        if(rows != cols){
            throw std::invalid_argument("Identity matrix must be a square!");
        }
        Matrix result(rows, cols);
        for(int i = 0; i < rows; ++i){
            for(int j = 0; j < cols; ++j){
                if(i == j){
                    result.data[i][j] = 1.0;
                }else{
                    result.data[i][j] = 0;
                }
            }
        }
        return result;
    }

    Matrix transpose() const{
        Matrix result(cols, rows);
        for(int i = 0; i < rows; ++i){
            for(int j = 0; j < cols; ++j){
                result.data[j][i] = data[i][j];
            }
        }
        return result;
    }

    Matrix inverse() const{
        if (rows != cols) throw std::invalid_argument("Matrix must be quadratic!");
        if(std::abs(determinant()) < 1e-12) throw std::runtime_error("Matrix is singular (det = 0)!");
        Matrix left = *this;
        Matrix right = left.identity_matrix();

        for(int i = 0; i < rows; ++i){
            // Search for pivot just like in gauss (what is the biggest value in the current column? => swap it)
            int pivot = i;
            for(int k = i + 1; k < rows; ++k){
                if(std::abs(left.data[k][i]) > std::abs(left.data[pivot][i])) pivot = k;
            }
            std::swap(left.data[i], left.data[pivot]);
            std::swap(right.data[i], right.data[pivot]); // swap in both matrices
            // Diagonal has to be the 1 (divide the whole row by the element on the diagonal)
            double diagVal = left.data[i][i];
            for(int j = 0; j < rows; ++j){
                left.data[i][j] /= diagVal;
                right.data[i][j] /= diagVal;
            }

            //eliminate up and down
            for(int k = 0; k < rows; ++k){
                if (k != i){ // everything to 0 except the pivot element
                    double factor = left.data[k][i];
                    for(int j = 0; j < rows; ++j){
                        left.data[k][j] -= factor * left.data[i][j];
                        right.data[k][j] -= factor * right.data[i][j];
                    }
                }
            }
        }
        return right;
    }

    Matrix pow(int a){
        if(rows != cols) throw std::invalid_argument("Calculation requires a square matrix! ");
        if(a < 0) throw std::invalid_argument("Negative exponents need the Inverse!");
        Matrix result(rows, cols);
        for(int i = 0; i < rows; ++i) result.data[i][i] = 1.0;

        for(int i = 0; i < a; ++i){
            result = result * (*this);
        }
        return result;
    }

    // help method for laplace expansion
    Matrix get_minor(int row_to_skip, int col_to_skip) const {
        Matrix submatrix(rows-1, cols-1);
        int sub_i = 0;
        for(int i = 0; i < rows; ++i){
            if(i == row_to_skip)continue;//skip this row

            int sub_j = 0;
            for(int j = 0; j < cols; ++j){
                if(j == col_to_skip) continue; //skip this column

                submatrix.data[sub_i][sub_j] = data[i][j];
                sub_j++;
            }
            sub_i++;
        }
        return submatrix;
    }

    //determinant
    double determinant() const {
        if(rows != cols){
            throw std::invalid_argument("Determinant requires a square matrix! ");
        }
        if(rows == 1){
            return data[0][0];
        }
        if(rows == 2){
            return(data[0][0] * data[1][1]) - (data[0][1] * data[1][0]);
        }
        if(rows == 3){
            return ((data[0][0] * data[1][1] * data[2][2]) + 
            (data[0][1] * data[1][2] * data[2][0]) + 
            (data[0][2] * data[1][0] * data[2][1])) - 
            ((data[2][0] * data[1][1] * data[0][2]) + 
            (data[2][1] * data[1][2] * data[0][0]) + 
            (data[2][2] * data[1][0] * data[0][1]));
        }
        // Laplace Expansion (4x4 to 6x6)
        // recursive approach: expanding along the first row
        if(rows >= 4 && rows <= 6){
            double det = 0;
            int sign = 1;

            for(int j = 0; j < cols; ++j){
                Matrix sub = get_minor(0, j); // get the minor matrix by removing row 0 and column j
                //recursive call to get determinant()
                det += sign * data[0][j] * sub.determinant();
                sign = -sign;
            }
            return det;
        } 
        // Gauss elimination to reduce time of calculating: O(n^3)
        int swapCount = 0;
        Matrix stepped = gauss(&swapCount); // transform to upper triangular matrix

        double det = 1.0;
        for (int i = 0; i < rows; ++i){
            det *= stepped.data[i][i];
        }
        if(swapCount % 2 != 0) det *= -1.0;
        // result very close to zero is treated as zero
        return (std::abs(det) < 1e-12) ? 0.0 : det;
    }

    Matrix gauss(int* swaps = nullptr) const{
        Matrix result = *this;
        if (swaps) *swaps = 0;
        for(int i = 0; i < rows; ++i){
            int pivotRow = i;
            // search for the biggest element in the current column => better numerical stability
            for(int k = i + 1; k < rows; ++k){
                if(std::abs(result.data[k][i]) > std::abs(result.data[pivotRow][i])){
                    pivotRow = k;
                }
            }
            if(pivotRow != i){
                // bigger element found => swap rows
                std::swap(result.data[i], result.data[pivotRow]);
                if(swaps) (*swaps)++;
            }
            // check if pivot element is not zero
            if(std::abs(result.data[i][i]) > 1e-15){
                for(int k = i + 1; k < rows; ++k){
                    //calculate the factor to create a zero in the current column
                    double factor = result.data[k][i] / result.data[i][i];
                    for(int j = i; j < rows; ++j){
                        // row - factor * row
                        result.data[k][j] -= factor * result.data[i][j];
                    }
                }
            }
        }
        return result;
    }



Matrix rref() const {
    Matrix res = *this;
    int r = 0;//pivot row index
    for (int c = 0; c < cols && r < rows; ++c) {//find pivots
        int pivotRow = r;
        for (int k = r + 1; k < rows; ++k) {//find largest value in this column
            if (std::abs(res.data[k][c]) > std::abs(res.data[pivotRow][c])) pivotRow = k;
        }
        //if column is nearly zero, skip to the next column
        if (std::abs(res.data[pivotRow][c]) < 1e-12) continue;
        //swap row with the row containing the largest pivot element
        std::swap(res.data[r], res.data[pivotRow]);
        //scale it to 1
        double diagVal = res.data[r][c];
        for (int j = c; j < cols; ++j) res.data[r][j] /= diagVal;
        //eliminate all elements (zero) above and below the current pivot
        for (int k = 0; k < rows; ++k) {
            if (k != r) {//skip the current pivot row itself
                double factor = res.data[k][c];
                for (int j = c; j < cols; ++j) res.data[k][j] -= factor * res.data[r][j];
            }
        }
        r++;
    }
    return res;
}

    vector<Matrix> basis() const{
        Matrix stepped = this->rref();
        vector<Matrix> basis;
        int currentRow = 0;

        for(int j = 0; j < cols && currentRow < rows; ++j){
            if(std::abs(stepped.data[currentRow][j]) > 1e-12){//search for pivot element
                Matrix columnVector(rows, 1);//create one vector that will be one of the basis vectors if we have a linearly independet vector (pivot after gauss)
                for(int i = 0; i < rows; ++i){
                    columnVector.data[i][0] = this->data[i][j];//fill the basis vector with the matrix from beginning
                }
                basis.push_back(columnVector);
                currentRow++;
            }
        }
        return basis;
    }

vector<Matrix> kernel() const {
    Matrix rref = this->rref();
    vector<Matrix> basis;
    
    // find out which column has a pivot
    // pivotRows[c] saves the row, where there is a pivot for column c (-1 if it is not a pivot)
    vector<int> pivotRows(cols, -1);
    int currentRow = 0;
    for (int c = 0; c < cols && currentRow < rows; ++c) {
        if (std::abs(rref.data[currentRow][c] - 1.0) < 1e-12) {
            pivotRows[c] = currentRow;
            currentRow++;
        }
    }

    // now just find out which are the empty columns (pivotRows[c] == -1)
    for (int c = 0; c < cols; ++c) {
        if (pivotRows[c] == -1) {
            Matrix v(cols, 1);
            v.data[c][0] = 1.0; // set free variables

            // calculate pivot variables (negative values of the rref)
            for (int p = 0; p < cols; ++p) {
                if (pivotRows[p] != -1) {
                    v.data[p][0] = -rref.data[pivotRows[p]][c];
                }
            }
            basis.push_back(v);
        }
    }
    return basis;
}

    void print() const {
        cout << "\n" << rows << "x" << cols << "):" << endl;
        for(int i = 0; i < rows; ++i){
            for(int j = 0; j < cols; ++j){
                if(std::abs(data[i][j]) < 1e-12){
                    cout << "0\t";
                }else{
                    cout << data[i][j] << "\t";
                }
            }
            cout << endl;
        }
        cout << endl;
    }

    void fill_user_matrix(){
        cout << "Enter your values for a " << rows << "x" << cols << " matrix: " << endl;
        for(int i = 0; i < rows; ++i){
            for(int j = 0; j < cols; ++j){
                cout << "element [" << i << "][" << j << "]: ";
                std::cin >> data[i][j];
            }
        }
        cout << endl;
    }


};

    Matrix create_matrix() {
        int r, c;
        cout << "How many rows? "; std::cin >> r;
        cout << "How many columns? "; std::cin >> c;

        Matrix m(r, c);
        m.fill_user_matrix();
        return m;

    }


int main(){
    int choice;
    cout << "-- calc for matrices --" << endl;
    cout << "1. Add matrices" << endl;
    cout << "2. Substract matrices" << endl;
    cout << "3. determinant" << endl;
    cout << "4. Gauss" << endl;
    cout << "5. Matrix multiplication (number/ 2nd matrix)" << endl;
    cout << "6. Matrix to the power" << endl;
    cout << "7. Transpose" << endl;
    cout << "8. Inverse" << endl;
    cout << "9. Rank, Basis and Kernel" << endl;
    cout << "End Program" << endl;
    cout << "Your choice? ";
    std::cin >> choice;

    switch (choice) {
        case 1: {
            cout << "\ndata for matrix a: " << endl;
            Matrix a = create_matrix();
            cout << "\ndata for matrix b: " << endl;
            Matrix b = create_matrix();

            Matrix res = a + b;
            res.print();
            break;
        }
        case 2: {
            cout << "\ndata for matrix a: " << endl;
            Matrix a = create_matrix();
            cout << "\ndata for matrix b: " << endl;
            Matrix b = create_matrix();

            Matrix res = a - b;
            res.print();
            break;

        }
        case 3: {
            cout << "\ndata for matrix: " << endl;
            Matrix a = create_matrix();
            cout << "determinant is: " << a.determinant() << endl;
            break;

        }
        case 4:{
            cout << "\ndata for matrix: " << endl;
            Matrix a = create_matrix();
            Matrix result = a.gauss();
            result.print();
            break;
        }
        case 5: {
            cout << "Multiplication with second matrix/vector (1) or one number (2)?" << endl;
            int v_choice;
            cout << "Your choice?" << endl;
            std::cin >> v_choice;
            switch (v_choice){
                case 1: { 
                    cout << "\ndata for matrix a: " << endl;
                    Matrix a = create_matrix();
                    cout << "\ndata for matrix b: " << endl;
                    Matrix b = create_matrix();
                    Matrix result = a * b;
                    result.print();
                    break;
                }
                case 2: {
                    cout << "\ndata for matrix a: " << endl;
                    Matrix a = create_matrix();
                    double number_m;
                    cout << "\nNumber to multiply matrix with?" << endl;
                    std::cin >> number_m;
                    Matrix result = a * number_m;
                    result.print();
                    break;
                }
            }
        }
        case 6:{
            cout << "\ndata for matrix: " << endl;
            Matrix a = create_matrix();
            cout << "What is the exponent?" << endl;
            int exp;
            std::cin >> exp;
            Matrix result = a.pow(exp);
            result.print();

        }
        case 7: {
            cout << "\ndata for matrix: " << endl;
            ((create_matrix()).transpose()).print();
            break;
        }
        case 8:{
            cout << "\ndata for matrix: " << endl;
            ((create_matrix()).inverse()).print();
            break;
        }
        case 9: {
            cout << "\ndata for matrix: " << endl;
            Matrix a = create_matrix();
            vector<Matrix> basis = a.basis();
            int rank = basis.size();
            cout << "Rank: " << rank << endl << endl;
            if(rank == 0){
                cout << "The basis is empty (Zero matrix)." << endl;
            } else{
                cout << "Basis: " << endl << endl;
                for(int i = 0; i < rank; ++i){
                    cout << "Vector " << (i + 1) << ":";
                    basis[i].print();// matrix printing method
                }
            }
            vector<Matrix> kernel = a.kernel();
            if(kernel.size() == 0){
                cout << "The matrix does not have any vectors that always result in zero => no kernel!" << endl;
            }else{
                cout << "Kernel is: " << endl << endl;
                for(size_t i = 0; i < kernel.size(); ++i){
                    cout << "Vector " << (i + 1) << endl;
                    kernel[i].print();
                }
            }
            
        }
        default:
            cout << "Invalid choice or program ended." << endl;
            break;
    }
}