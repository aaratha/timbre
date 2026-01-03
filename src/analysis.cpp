#include <iostream>

#include "analysis.hpp"
#include "globals.hpp"
#include "miniaudio.h"

AnalysisCore::AnalysisCore() {}

AnalysisCore::~AnalysisCore() {}

void AnalysisCore::readFile(const std::string &filename) {
  ma_decoder decoder;

  ma_decoder_config config =
      ma_decoder_config_init(ma_format_f32, 2, DEVICE_SAMPLE_RATE);

  // 0 = use file's channels and sample rate
  ma_result result = ma_decoder_init_file(filename.c_str(), &config, &decoder);
  if (result != MA_SUCCESS) {
    std::cerr << "Failed to open file\n";
    return;
  }

  // Get total frame count
  ma_uint64 totalFrameCount = 0;
  result = ma_decoder_get_length_in_pcm_frames(&decoder, &totalFrameCount);
  if (result != MA_SUCCESS) {
    std::cerr << "Failed to read file length\n";
    ma_decoder_uninit(&decoder);
    return;
  }

  // Prepare buffer: channels * frames
  std::vector<float> buffer(totalFrameCount * decoder.outputChannels);

  // Read frames into buffer
  ma_uint64 framesRead = 0;
  result = ma_decoder_read_pcm_frames(&decoder, buffer.data(), totalFrameCount,
                                      &framesRead);
  if (result != MA_SUCCESS || framesRead != totalFrameCount) {
    std::cerr << "Warning: did not read all frames\n";
  }

  // store input sample raw data
  inputRaw = buffer;

  std::cout << "Read " << framesRead << " frames from " << filename << "\n";

  ma_decoder_uninit(&decoder);
}

std::vector<float> &AnalysisCore::getInputRaw() { return inputRaw; }
