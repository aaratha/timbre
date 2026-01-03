#include "audio.hpp"
#include "globals.hpp"

#include <iostream>

AudioCore::AudioCore(AnalysisCore &analysisCore)
    : analysisCore(analysisCore) {
  if (audioInitialized)
    return;

  deviceConfig = ma_device_config_init(ma_device_type_playback);
  deviceConfig.playback.format = DEVICE_FORMAT;
  deviceConfig.playback.channels = DEVICE_CHANNELS;
  deviceConfig.sampleRate = static_cast<ma_uint32>(DEVICE_SAMPLE_RATE);
  deviceConfig.dataCallback = AudioCore::dataCallback;
  deviceConfig.pUserData = this;

  osc.sampleRate = DEVICE_SAMPLE_RATE;

  if (ma_device_init(NULL, &deviceConfig, &device) != MA_SUCCESS) {
    std::cerr << "ma_device_init failed\n";
    return;
  }

  if (ma_device_start(&device) != MA_SUCCESS) {
    std::cerr << "ma_device_start failed\n";
    ma_device_uninit(&device);
    return;
  }

  playhead.store(0);
  audioInitialized = true;
  running.store(true);
}

AudioCore::~AudioCore() {
  if (audioInitialized) {
    ma_device_uninit(&device);
    audioInitialized = false;
    running.store(false);
  }
}

void AudioCore::dataCallback(ma_device *pDevice, void *pOutput,
                             const void * /*pInput*/, ma_uint32 frameCount) {
  auto *core = static_cast<AudioCore *>(pDevice->pUserData);
  float *out = static_cast<float *>(pOutput);

  if (!core->callbackSeen.exchange(true)) {
    std::cerr << "audio callback active\n";
  }

  core->processAudio(out, frameCount);
}


void AudioCore::processAudio(float *out, ma_uint32 frameCount) {
  for (ma_uint32 frame = 0; frame < frameCount; ++frame) {
    if (!analysisCore.getInputRaw().empty()) {
      // use playhead to track position in inputRaw
      size_t idx = playhead.load();
      size_t channels = static_cast<size_t>(DEVICE_CHANNELS);
      size_t total = analysisCore.getInputRaw().size();
      for (ma_uint32 ch = 0; ch < DEVICE_CHANNELS; ++ch) {
        size_t sampleIndex = (idx + ch) % total;
        *out++ = analysisCore.getInputRaw()[sampleIndex];
      }
      idx = (idx + channels) % total;
      playhead.store(idx);
    } else {
      // Generate a simple sine wave if no input file is loaded
      float sample = 0.2f * osc.process();
      for (ma_uint32 ch = 0; ch < DEVICE_CHANNELS; ++ch) {
        *out++ = sample;
      }
    }
  }
}
