#include <iostream>
#include "ParseException.h"
#include "Scene.h"

using namespace glm;

Settings parseArgs(int argc, char **argv) {
    int paramCount = 4;
    if (argc == 1)
        return Settings();
    if (argc < 2 * paramCount - 1) {
        throw ParseException(ParseException::NOT_ENOUGH);
    }
    Settings settings("", 0, 0);
    for (int i = 1; i < 2 * paramCount - 1; i += 2) {
        std::string arg = argv[i];
        if (arg == "-out") {
            if (settings.out.length() > 0) {
                throw ParseException(ParseException::SAME_ARGS);
            }
            settings.out = argv[i + 1];
        } else if (arg == "-scene") {
            if (settings.scene != 0) {
                throw ParseException(ParseException::SAME_ARGS);
            }
            settings.scene = atoi(argv[i + 1]);
        } else if (arg == "-threads") {
            if (settings.threads != 0) {
                throw ParseException(ParseException::SAME_ARGS);
            }
            settings.threads = atoi(argv[i + 1]);
        }
    }
    return settings;
}

int main(int argc, char **argv) {
    Settings settings;
    try {
        settings = parseArgs(argc, argv);
    } catch (ParseException &e) {
        e.printMessage();
        return -1;
    }
    if (settings.scene != 1)
        return 0;
    Scene scene;
    return scene.run(settings);
}