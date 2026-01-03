#include <QGuiApplication>
#include <QQmlApplicationEngine>

#include "audio.hpp"
#include "ui.hpp"
#include "timbre.hpp"

int main(int argc, char *argv[]) {
    QGuiApplication app(argc, argv);

    QQmlApplicationEngine engine;

    AnalysisCore analysisCore;
    analysisCore.readFile(argv[1]); // Load audio file from command line argument
    
    AudioCore audioCore(analysisCore);
    
    engine.loadFromModule("MyApp", "MainView");

    if (engine.rootObjects().isEmpty())
        return -1;

    return app.exec();
}
