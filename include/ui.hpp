#pragma once

#include <QObject>
#include <QVariant>

class AnalysisCore;
class AudioCore;

class UiController : public QObject {
  Q_OBJECT
  Q_PROPERTY(int lastClicked READ lastClicked NOTIFY lastClickedChanged)
  Q_PROPERTY(int binCount READ binCount NOTIFY binCountChanged)
  Q_PROPERTY(QVariantList umapPoints READ umapPoints NOTIFY umapPointsChanged)

public:
  explicit UiController(AudioCore *audioCore, AnalysisCore *analysisCore,
                        QObject *parent = nullptr);

  Q_INVOKABLE void rectangleClicked(int index);
  int lastClicked() const;
  int binCount() const;
  QVariantList umapPoints() const;

signals:
  void lastClickedChanged();
  void binCountChanged();
  void umapPointsChanged();

private:
  AudioCore *audioCore;
  AnalysisCore *analysisCore;
  int lastClickedIndex{-1};
};
