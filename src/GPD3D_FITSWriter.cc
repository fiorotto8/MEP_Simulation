// GPD3D_FITSWriter.cc

#include "GPD3D_FITSWriter.hh"

void GPD3D_FITSWriter::writeFITS(const char* fileName,
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
                                  const std::vector<std::vector<std::vector<double>>>& Maptoa3D) {
                                  // const std::vector<std::vector<std::vector<double>>>& Maptot3D) {
    fitsfile *fptr;
    int status = 0;

    // Check if the file exists and delete it
    if (std::remove(fileName) != 0) {
        // Error handling if the file deletion fails (optional)
        std::perror("Error deleting existing file");
    }

    // Create a new file
    fits_create_file(&fptr, fileName, &status);

    // Check if there was an error creating the file
    if (status) {
        fits_report_error(stderr, status);
        return;
    }

    const char *ttype[] = {"EventID", "PhoEnergy", "AbsPos_X", "AbsPos_Y", "AbsPos_Z",
                           "PEEnergy", "PE_PHI", "PE_THETA", "AUGEnergy", "AUG_PHI", "AUG_THETA", "TrkLength", "Ion_Pos_X", "Ion_Pos_Y", "Ion_Pos_Z"};
    const char *tform[] = {"J", "D", "D", "D", "D", "D", "D", "D", "D", "D", "D", "D", "PD", "PD", "PD"};
    const char *tunit[] = {"", "keV", "mm", "mm", "mm", "keV", "deg", "deg", "keV", "deg", "", "mm", "", "", ""};
    fits_create_tbl(fptr, BINARY_TBL, 0, 15, (char**)ttype, (char**)tform, (char**)tunit, "MonteCarlo", &status);

    long fpixel = 1;
    long naxis2 = eventIDList.size();
    
    fits_write_col(fptr, TINT, 1, fpixel, 1, naxis2, const_cast<int*>(eventIDList.data()), &status);
    fits_write_col(fptr, TDOUBLE, 2, fpixel, 1, naxis2, const_cast<double*>(phoEngList.data()), &status);
    fits_write_col(fptr, TDOUBLE, 3, fpixel, 1, naxis2, const_cast<double*>(absPosXList.data()), &status);
    fits_write_col(fptr, TDOUBLE, 4, fpixel, 1, naxis2, const_cast<double*>(absPosYList.data()), &status);
    fits_write_col(fptr, TDOUBLE, 5, fpixel, 1, naxis2, const_cast<double*>(absPosZList.data()), &status);
    fits_write_col(fptr, TDOUBLE, 6, fpixel, 1, naxis2, const_cast<double*>(peEngList.data()), &status);
    fits_write_col(fptr, TDOUBLE, 7, fpixel, 1, naxis2, const_cast<double*>(pePhiList.data()), &status);
    fits_write_col(fptr, TDOUBLE, 8, fpixel, 1, naxis2, const_cast<double*>(peThetaList.data()), &status);
    fits_write_col(fptr, TDOUBLE, 9, fpixel, 1, naxis2, const_cast<double*>(augEngList.data()), &status);
    fits_write_col(fptr, TDOUBLE, 10, fpixel, 1, naxis2, const_cast<double*>(augPhiList.data()), &status);
    fits_write_col(fptr, TDOUBLE, 11, fpixel, 1, naxis2, const_cast<double*>(augThetaList.data()), &status);
    fits_write_col(fptr, TDOUBLE, 12, fpixel, 1, naxis2, const_cast<double*>(trkLenList.data()), &status);

    // Write Ion Pos X, Y, Z
    for (size_t i = 0; i < ionPosXList.size(); ++i) {
        fits_write_col(fptr, TDOUBLE, 13, fpixel + i, 1, ionPosXList[i].size(), const_cast<double*>(ionPosXList[i].data()), &status);
        fits_write_col(fptr, TDOUBLE, 14, fpixel + i, 1, ionPosYList[i].size(), const_cast<double*>(ionPosYList[i].data()), &status);
        fits_write_col(fptr, TDOUBLE, 15, fpixel + i, 1, ionPosZList[i].size(), const_cast<double*>(ionPosZList[i].data()), &status);
    }

    // Flatten gridResults3D and write it to a new table in the existing HDU
    // Create a new table for gridResults3D
    // const char *ttype_grid[] = {"DIFFUSE_POS_X", "DIFFUSE_POS_Y", "DIFFUSE_POS_Z", "COL_LIST", "ROW_LIST", "COUNTING_MAP", "MAP_TOA", "MAP_TOT"};
    const char *ttype_grid[] = {"COL_LIST", "ROW_LIST", "COUNTING_MAP", "MAP_TOA"};
    // const char *tform_grid[] = {"PD", "PD", "PD", "PD", "PD",  "PI", "PD", "PD"};
    const char *tform_grid[] = {"PD", "PD",  "PD", "PD"};
    // const char *tunit_grid[] = {"", "", "", "", "", "", "", ""};
    const char *tunit_grid[] = {"", "", "", ""};
    fits_create_tbl(fptr, BINARY_TBL, 0, 4, const_cast<char**>(ttype_grid), const_cast<char**>(tform_grid), const_cast<char**>(tunit_grid), "DATA", &status);

    // Write Gauss Pos X, Y, Z
    // for (size_t i = 0; i < gaussPosXList.size(); ++i) {
    //     fits_write_col(fptr, TDOUBLE, 1, fpixel + i, 1, gaussPosXList[i].size(), const_cast<double*>(gaussPosXList[i].data()), &status);
    //     fits_write_col(fptr, TDOUBLE, 2, fpixel + i, 1, gaussPosYList[i].size(), const_cast<double*>(gaussPosYList[i].data()), &status);
    //     fits_write_col(fptr, TDOUBLE, 3, fpixel + i, 1, gaussPosZList[i].size(), const_cast<double*>(gaussPosZList[i].data()), &status);
    // }

    for (size_t i = 0; i < mapcolList.size(); ++i) {
        fits_write_col(fptr, TDOUBLE, 1, fpixel + i, 1, mapcolList[i].size(), const_cast<double*>(mapcolList[i].data()), &status);
        fits_write_col(fptr, TDOUBLE, 2, fpixel + i, 1, maprowList[i].size(), const_cast<double*>(maprowList[i].data()), &status);
    }

    for (size_t i = 0; i < countingMap3D.size(); ++i) {
        std::vector<int> flatCountingMap;
        for (const auto& row : countingMap3D[i]) {
            flatCountingMap.insert(flatCountingMap.end(), row.begin(), row.end());
        }

        // Write the flattened gridResults3D to the new table
        fits_write_col(fptr, TINT, 3, fpixel + i, 1, flatCountingMap.size(), const_cast<int*>(flatCountingMap.data()), &status);
    }

    for (size_t i = 0; i < Maptoa3D.size(); ++i) {
        std::vector<double> flatMapToA;
        for (const auto& row : Maptoa3D[i]) {
            flatMapToA.insert(flatMapToA.end(), row.begin(), row.end());
        }
        // Write the flattened gridResults3D to the new table
        fits_write_col(fptr, TDOUBLE, 4, fpixel + i, 1, flatMapToA.size(), const_cast<double*>(flatMapToA.data()), &status);
    }

    // for (size_t i = 0; i < Maptot3D.size(); ++i) {
    //     std::vector<double> flatMapToT;
    //     for (const auto& row : Maptot3D[i]) {
    //         flatMapToT.insert(flatMapToT.end(), row.begin(), row.end());
    //     }

    //     // Write the flattened gridResults3D to the new table
    //     fits_write_col(fptr, TDOUBLE, 8, fpixel + i, 1, flatMapToT.size(), const_cast<double*>(flatMapToT.data()), &status);
    // }
    
    fits_close_file(fptr, &status);

    if (status) {
        fits_report_error(stderr, status);
        return;
    }
}
