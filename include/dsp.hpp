#pragma once

#include <vector>

#include "audio.hpp"

namespace DSP {
// Compute power spectrum from input sample
void computePS(const float *input, std::vector<std::complex<float>> &output,
               int size);

// Compute MFCCs from power spectrum
void computeMFCC(const float *input, std::vector<float> &output, int size);

// Compute centroid from input sample
void computeCentroid(const float *input, int size, float &centroid);

// Compute flux from input sample
void computeFlux(const float *input, const float *prevInput, int size,
                 float &flux);

// Compute rolloff from input sample
void computeRolloff(const float *input, int size, float &rolloff);

// Compute UMAP from MFCCs, centroid, flux, and rolloff
void computeUMAP(const std::vector<std::vector<float>> &input,
                 std::vector<std::vector<float>> &output, int nComponents);
}; // namespace DSP
