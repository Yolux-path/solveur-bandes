// Fourni par l'enseignant (un #include ajouté).
#include <string>
#include <iostream> // The std library for i/o
#include <fstream>  // For file streams
#include <sstream>  // For string streams
#include <iomanip>  // Manipulate i/o, e.g. formatting
#include <cmath>
#include <cstdlib>

#include "mat_vec.h" // My struct for matrices

// Read the coefficients of the matrix from a text file
matrix readMatrix(std::string fileName) {
    matrix M;
    // We setup an input stream from the file
    std::ifstream mFile(fileName);
    
    if (not mFile.good()) {
        std::cerr << "Unable to open file " << fileName << std::endl;
        abort();
    }
    
    // We reach the first line that does not start with '%' or '#'
    std::string line;
    getline(mFile,line);
    while (line.front() == '%' || line.front() == '#')
        getline(mFile,line);
    
    // Now, line contains the numbers of row and columns of the matrix
    std::istringstream iss(line);
    iss >> M.n >> M.m;
    if (M.n!=M.m) {
        std::cerr << "Error: this program works only with square matrices."
                  << std::endl;
        abort();
    }
        
    // Set the size of the matrix -- allocate memory
    M.a.resize(M.n);
    for (uint i=0u; i<M.n; i++) M.a[i].resize(M.m);

    // Now the file stream is located at the first row of coefficients
    for (uint i=0u; i<M.n; i++) {
        for (uint j=0u; j<M.m; j++) mFile >> M.a[i][j];
    }
    
    return M;
}

// Read the coefficients of the matrix from a text file
mat_band readMatrixBand(std::string fileName) {
    mat_band M;
    // We setup an input stream from the file
    std::ifstream mFile(fileName);
    
    if (not mFile.good()) {
        std::cerr << "Unable to open file " << fileName << std::endl;
        abort();
    }
    
    // We reach the first line that does not start with '%' or '#'
    std::string line;
    getline(mFile,line);
    while (line.front() == '%' || line.front() == '#')
        getline(mFile,line);
    
    // Now, line contains the numbers of row and columns of the matrix
    std::istringstream iss(line);
    iss >> M.n;
    // Next line : kl and ku
    mFile >> M.kl >> M.ku;
        
    // Set the size of the matrix -- allocate memory
    M.col.resize(M.n);
    for (uint j=0u; j<M.n; j++) M.col[j].resize(1u+M.kl+M.ku);
    
    // Now the file stream is located at the first row of coefficients
    for (uint k=0u; k<1u+M.kl+M.ku; k++)
        for (uint j=0u; j<M.n; j++) mFile >> M.col[j][k];
        
    return M;
}

// Reads coefficients of a sparse matrix from a text file following the .mtx
// format from MatrixMarket:
// https://math.nist.gov/MatrixMarket/formats.html#MMformat
// and store it as a band matrix
mat_band readMatrixBandFromMtx(std::string fileName) {
    mat_band M;
    // We setup an input stream from the file
    std::ifstream mFile(fileName);
    
    if (not mFile.good()) {
        std::cerr << "Unable to open file " << fileName << std::endl;
        abort();
    }
    
    // We reach the first line that does not start with '%' or '#'
    std::string line;
    getline(mFile,line);
    while (line.front() == '%' || line.front() == '#')
        getline(mFile,line);
    
    // Now line contains rows, columns, and nnz
    uint n,p,nnz;
    std::istringstream iss(line);
    iss >> n >> p >> nnz;

    if (n!=p) {
        std::cerr
            << "Error: this program works only with square matrices."
            << std::endl;
        abort();
    }
        
    // Set the size of the tmp vectors 
    std::vector<uint> rows(nnz);
    std::vector<uint> cols(nnz);
    std::vector<double> vals(nnz);
    
    // Now the file stream is located at the definitions of the coefficients
    // Look for kl and ku in the same loop
    M.kl = 0;
    M.ku = 0;  
    for (uint k=0u; k<nnz; k++) {
        mFile >> rows[k] >> cols[k] >> vals[k];
        // Indices of MatrixMarket start at 1
        rows[k]--;
        cols[k]--;
        M.ku = std::max(int(M.ku),int(cols[k])-int(rows[k]));
        M.kl = std::max(int(M.kl),-(int(cols[k])-int(rows[k])));
    }

    // We have to convert these arrays to a mat_band storage.
    M.n = n;
    M.col.resize(n);
    for (uint j=0u; j<M.n; j++) M.col[j].resize(1u+M.kl+M.ku, 0.e0);
    // Add values
    for (uint i=0u; i<nnz; i++) {
        uint k = M.ku - (int(cols[i]) - int(rows[i])); // k := ku - (j-i)
        M.col[cols[i]][k] += vals[i];
    }
    
    return M;
}

// Reads the coefficients of a vector from a file
vec readVector(std::string fileName) {
    vec B;
    // We setup an input stream from the file
    std::ifstream mFile(fileName);
    
    if (not mFile.good()) {
        std::cerr << "Unable to open file " << fileName << std::endl;
        abort();
    }
    
    // We reach the first line that does not start with '%' or '#'
    std::string line;
    getline(mFile,line);
    while (line.front() == '%' || line.front() == '#')
        getline(mFile,line);
    
    // Now, line contains the size of the vector
    uint n;
    std::istringstream iss(line);
    iss >> n;
        
    // Set the size of the matrix  
    B.resize(n);

    // Now the file stream is located at the first coefficient
    for (uint i=0u; i<n; i++) mFile >> B[i];

    return B;
}

// Reads coefficients of a vector (RHS) from a text file following the .mtx
// format from MatrixMarket:
// https://math.nist.gov/MatrixMarket/formats.html#MMformat
vec readVectorFromMtx(std::string fileName) {

    // We setup an input stream from the file
    std::ifstream mFile(fileName);
    
    if (not mFile.good()) {
        std::cerr << "Unable to open file " << fileName << std::endl;
        abort();
    }
    
    // We reach the first line that does not start with '%' or '#'
    std::string line;
    getline(mFile,line);
    while (line.front() == '%' || line.front() == '#')
        getline(mFile,line);
    
    // Now line contains rows, columns
    uint n,p;
    std::istringstream iss(line);
    iss >> n >> p;

    if (p!=1) {
        std::cerr
            << "Error: this program works with only single RHS vectors."
            << std::endl;
        abort();
    }

    // Set the size of the vector
    vec X(n,0.0e0);
    
    // Now the file stream is located at the definitions of the coefficients
    for (uint k=0u; k<n; k++) mFile >> X[k];
    
    return X;
}

// Prints on screen the contents of a matrix
void displayMatrix(const matrix& M) {

    std::cout << "Matrix of size: " << M.n << ", " << M.m << std::endl;
    for (uint i=0u; i<M.n; i++) {
        for (uint j=0u; j<M.m; j++)
            std::cout << std::setw(14) << std::scientific << M.a[i][j];
        std::cout << std::endl;        
    }
}

// Prints on screen the contents of a matrix
void displayMatrixBand(const mat_band& M) {

    double a;
    int k;
    
    std::cout << "Band matrix of size " << M.n
              << " with bandwidth kl, ku = "
              << M.kl << " , " << M.ku << std::endl;
    
    for (int i=0; i<int(M.n); i++) {
        for (int j=0; j<int(M.n); j++) {
            // note: k = ku - (j-i) is the row of entry (i,j)
            k = M.ku - (j-i);
            if ((j-i)>int(M.ku) or (j-i)<-int(M.kl))
                a = 0.;
            else
                a = M.col[j][k];
            std::cout << std::scientific << std::setw(14) << a << " ";
        }
        std::cout << std::endl;        
    }
}

// Write matrix to an ofstream
void writeMatrix(const matrix& M, std::ofstream& output_file) {

    output_file << "Matrix of size: " << M.n << ", " << M.m << std::endl;
    for (uint i=0u; i<M.n; i++) {
        for (uint j=0u; j<M.m; j++)
            output_file << std::setw(14) << std::scientific << M.a[i][j];
        output_file << std::endl;        
    }
}

// Write the Band matrix to an ofstream
void writeMatrixBand(const mat_band& M, std::ofstream& output_file) {

    output_file << "# This is a square band matrix" << std::endl;
    output_file << M.n << std::endl; // size
    output_file << M.kl << " " << M.ku << std::endl; // kl and ku    

    for (uint k=0u; k<M.kl+1+M.ku; k++) {
        for (uint j=0u; j<M.n; j++)
            output_file << std::scientific << std::setw(14) << M.col[j][k];
        output_file << std::endl;
    }
}

// Prints on screen the contents of a vector
void displayVector(const vec& B) {

    uint n = B.size();
    
    std::cout << "Vector of size: " << n << std::endl;
    for (uint i=0u; i<n; i++) 
        std::cout << std::scientific << std::setw(14) << B[i];
    std::cout << std::endl;        
}

// Write the contents of a vector to an ofstream
void writeVector(const vec& B, std::ofstream& output_file) {

    uint n = B.size();
    
    output_file << "Vector of size: " << n << std::endl;
    for (uint i=0u; i<n; i++) output_file << B[i] << " ";
    output_file << std::endl;        
}

// Plot Phi and grad(Phi) = vector of 2D data on a cartesian grid to a
// txt file for gnuplot figures.
void vectorToGnuplot(const vec& Phi, const vec& vx, const vec& vy,std::string fileName, const unsigned int& sample) {

    unsigned int n2 = Phi.size();
    unsigned int n = sqrt(n2);

    double h = 1./(n-1u);
    double x,y;
    
    std::ofstream gp (fileName);

    gp << "# data on a nxn cartesian grid of (0,1)x(0,1)" << std::endl;
    for (unsigned int i=0; i<n; i+=sample) {
        x = i*h;
        for (unsigned int j=0; j<n; j+=sample) {
            y = j*h;
            gp << std::scientific << std::setw(14) << x << " "
               << std::setw(14) << y << " "
               << std::setw(14) << Phi[i+n*j]
               << std::setw(14) << vx[i+n*j]
               << std::setw(14) << vy[i+n*j]
               << std::endl;
        }
        gp << std::endl;
    }
}


