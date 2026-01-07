#include "ui.hpp"

#include <algorithm>
#include <iostream>

#include "analysis.hpp"
#include "audio.hpp"

UiController::UiController(AudioCore *audioCore, AnalysisCore *analysisCore,
                           QObject *parent)
    : QObject(parent), audioCore(audioCore), analysisCore(analysisCore) {
  emit binCountChanged();
  emit umapPointsChanged();
}

void UiController::rectangleClicked(int index) {
  lastClickedIndex = index;
  std::cerr << "rectangle clicked: " << index << "\n";
  emit lastClickedChanged();

  if (audioCore) {
    audioCore->setBinIndex(static_cast<size_t>(index));
  }
}

int UiController::lastClicked() const { return lastClickedIndex; }

int UiController::binCount() const {
  if (!analysisCore) {
    return 0;
  }
  return static_cast<int>(analysisCore->getBinCount());
}

QVariantList UiController::umapPoints() const {
  QVariantList points;
  if (!analysisCore) {
    return points;
  }

  const auto &xs = analysisCore->getBinTimbreX();
  const auto &ys = analysisCore->getBinTimbreY();
  const size_t count = std::min(xs.size(), ys.size());
  if (count == 0) {
    return points;
  }

  float minX = xs[0];
  float maxX = xs[0];
  float minY = ys[0];
  float maxY = ys[0];
  for (size_t i = 1; i < count; ++i) {
    minX = std::min(minX, xs[i]);
    maxX = std::max(maxX, xs[i]);
    minY = std::min(minY, ys[i]);
    maxY = std::max(maxY, ys[i]);
  }

  const float rangeX = maxX - minX;
  const float rangeY = maxY - minY;

  points.reserve(static_cast<int>(count));
  for (size_t i = 0; i < count; ++i) {
    const float normX = rangeX > 0.0f ? (xs[i] - minX) / rangeX : 0.5f;
    const float normY = rangeY > 0.0f ? (ys[i] - minY) / rangeY : 0.5f;
    QVariantMap point;
    point.insert("x", normX);
    point.insert("y", normY);
    point.insert("index", static_cast<int>(i));
    points.append(point);
  }

  return points;
}
