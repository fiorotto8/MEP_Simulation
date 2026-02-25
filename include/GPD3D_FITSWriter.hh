#ifndef GPD3D_FITSWRITER_HH
#define GPD3D_FITSWRITER_HH

#include <cstdio>
#include <vector>
#include <iostream>
#include "fitsio.h"

class GPD3D_FITSWriter {
public:
    void static writeFITS(const char* fileName,
                                  const std::vector<int>& eventIDList,
                                  const std::vector<double>& phoEngList,
                                  const std::vector<double>& absPosXList,
                                  const std::vector<double>& absPosYList,
                                  const std::vector<double>& absPosZList,
                                  const std::vector<double>& peEngList,
                                  const std::vector<double>& pePhiList,
                                  const std::vector<double>& peThetaList,
                                  const std::vector<double>& augEngList,
                                  const std::vector<double>& augPhiList,
                                  const std::vector<double>& augThetaList,
                                  const std::vector<double>& trkLenList,
                                  const std::vector<std::vector<double>>& ionPosXList,
                                  const std::vector<std::vector<double>>& ionPosYList,
                                  const std::vector<std::vector<double>>& ionPosZList,
                                  // const std::vector<std::vector<double>>& gaussPosXList,
                                  // const std::vector<std::vector<double>>& gaussPosYList,
                                  // const std::vector<std::vector<double>>& gaussPosZList,
                                  const std::vector<std::vector<double>>& mapcolList,
                                  const std::vector<std::vector<double>>& maprowList,
                                  const std::vector<std::vector<std::vector<int>>>& countingMap3D,
                                  const std::vector<std::vector<std::vector<double>>>& Maptoa3D);
                                  // const std::vector<std::vector<std::vector<double>>>& Maptot3D);
};

#endif // GPD3D_FITSWRITER_HH
