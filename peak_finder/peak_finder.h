#pragma once

#ifndef PEAKFINDER_H
#define PEAKFINDER_H

#include <vector>       // for vector
#include <algorithm>    // for min_element|max_element
#include "../../PlasmaDataProcessing_DLL/PlasmaDataProcessing/Matrix_Class/mynamespace.h"
#include "../error_codes.h"

using namespace myspace;

/* find peaks */
namespace PeakFinder {
    const double EPS = 2.2204E-16f;

    /*
        Inputs
        x0: input signal
        extrema: 1 if maxima are desired, -1 if minima are desired
        includeEndpoints - If true the endpoints will be included as possible extrema otherwise they will not be included
        Output
        peakInds: Indices of peaks in x0
    */
    int findPeaks(_In_ const std::vector <double>&,
                  _Out_ std::vector <int>&,
                  _In_opt_ const bool = false,
                  _In_opt_ const int = 1);
};

#endif