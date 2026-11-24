#include <fstream>
#include <iostream>
#include <string>

#include "mat_vec.h"      // Matrices (full and band) and vectors
#include "mat_vec_io.h"   // I/O for matrices and vectors
#include "mat_vec_norm.h" // Norms of matrices and vectors
#include "linalg_band.h"  // Linear algebra algorithms for band matrices
#include "table.h"        // Table for the report

int main(int argc, char* argv[]) {
    
    //Example of usage:
    //./mini_prog ./.../NS_4.txt ./.../result_NS_4.txt

    // Check the number of command-line arguments
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input_file> <output_file>" << std::endl;
        return EXIT_FAILURE;
    }


    // Get input and output file names from command-line arguments
    std::string inputFile = argv[1];
    std::string outputFile = argv[2];

    // Open the output file
    std::ofstream results(outputFile);
    if (!results.is_open()) {
        std::cerr << "Error opening the output file: " << outputFile << std::endl;
        return EXIT_FAILURE;
    }

    // Read a band matrix from the input file in Matrix Market format

    mat_band A = inputFile.ends_with(".mtx") ? readMatrixBandFromMtx(inputFile)
                                             : readMatrixBand(inputFile);

    // Write the results to the output file using functions from the 'table' module
    writeMatrixBandAtable(results, ComputeAInfos(A));
    writeMatrixBandBtable(results, ComputeBInfos(A));
    writeMatrixBandCtable(results, ComputeCInfos(A));
    writeMatrixBandDtable(results, ComputeDInfos(A));
    writeMatrixBandEtable(results, ComputeEInfos(A));
    writeMatrixBandFtable(results, ComputeFInfos(A));

    // Close the output file
    results.close();

    return EXIT_SUCCESS;
}
