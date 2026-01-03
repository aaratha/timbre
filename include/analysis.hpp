#pragma once

#include <vector>

class AnalysisCore {
  std::vector<float> inputRaw;
  size_t inputBinSize;
  std::vector<std::vector<float>> binnedInput;

public:
  AnalysisCore();
  ~AnalysisCore();

  void readFile(const std::string &filename);
  std::vector<float> &getInputRaw();
  void binInput(int binSize);
};
