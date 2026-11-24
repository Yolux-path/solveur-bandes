#ifndef __TABLE_INCLUDED__
#define __TABLE_INCLUDED__

#include "mat_vec.h"

typedef unsigned int uint;

std::vector<uint> ComputeAInfos(const mat_band& M);

std::vector<uint> ComputeBInfos(const mat_band& M);

std::vector<double> ComputeCInfos(const mat_band& M);

std::vector<double> ComputeDInfos(const mat_band& M);

std::vector<double> ComputeEInfos(const mat_band& M);

std::vector<double> ComputeFInfos(const mat_band& M);

void writeMatrixBandAtable(std::ofstream& output_file, const std::vector<uint>& Infos);

void writeMatrixBandBtable(std::ofstream& output_file, const std::vector<uint>& Infos);

void writeMatrixBandCtable(std::ofstream& output_file, const std::vector<double>& Infos);

void writeMatrixBandDtable(std::ofstream& output_file, const std::vector<double>& Infos);

void writeMatrixBandEtable(std::ofstream& output_file, const std::vector<double>& Infos);

void writeMatrixBandFtable(std::ofstream& output_file, const std::vector<double>& Infos);

#endif
