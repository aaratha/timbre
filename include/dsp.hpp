#pragma once

#include <complex>
#include <vector>

namespace DSP {

// Apply Hann window to input samples
void hannWindow(std::vector<float> &input, int size);

// Compute FFT of input samples (re + im)
void computeFFT(const std::vector<float> &input,
                std::vector<std::complex<float>> &output, int size);

// Compute STFT of input samples with given size and hop size
// Returns vector of frames, only unique frequency bins (size/2 + 1) for each
// frame (for analysis)
void computeSTFT(const std::vector<float> &input,
                 std::vector<std::vector<std::complex<float>>> &output,
                 int fftSize, int frameCount);

void computeIFFT(const std::vector<std::complex<float>> &input,
                 std::vector<float> &output, int size);

// Compute power spectrum from FFT input
void computePowerSpectrum(const std::vector<std::complex<float>> &input,
                          std::vector<float> &output, int size);

// input: Power spectrum (length = nSpec)
// output: Mel-scale spectral envelope (length = nMelBands) - NOT LOG-SCALED YET
void computeSpectralEnv(const std::vector<float> &powerSpec,
                        std::vector<float> &melEnv, int sampleRate, int fftSize,
                        int nMelBands);

// Compute MFCCs from power spectrum
// input: Mel-scale spectral envelope (length = nMelBands)
void computeMFCC(const std::vector<float> &input, std::vector<float> &output,
                 int size);

// Compute centroid from input power spectrum
void computeCentroid(const std::vector<float> &input, int size, int sampleRate,
                     float &centroid);

// Compute flux from current and previous power spectra
void computeFlux(const std::vector<float> &psCurr,
                 const std::vector<float> &psPrev, int size, float &flux);

// Compute rolloff from input power spectrum
// ps: power or magnitude spectrum (length = size/2 + 1)
// size: original FFT size
// sampleRate: audio sample rate
// threshold: fraction of total spectral energy to consider (0.85 = 85%)
// rolloff: output frequency in Hz
void computeRolloff(const std::vector<float> &ps, int size, int sampleRate,
                    float threshold, float &rolloff);

// Compute UMAP from MFCCs, centroid, flux, and rolloff
// Returns 2D coordinates for each input feature vector
void computeUMAP(const std::vector<std::vector<float>> &inputFeatures,
                 std::vector<float> &outputX, std::vector<float> &outputY);
}; // namespace DSP
