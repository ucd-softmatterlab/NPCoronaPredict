#ifndef CONFIG_FILE_READER__H__
#define CONFIG_FILE_READER__H__

#include "Config.h"
#include "StringFormat.h"

#include <string>
#include <vector>
#include <fstream>
#include <iostream>

class ConfigFileReader {
private:
    std::vector<std::string> m_switches;
    std::vector<std::string> m_keys;
    std::vector<std::string> m_values;

public:
    ConfigFileReader(const std::string& filename) {
        ReadConfigFile(filename);
    }

    void ReadConfigFile(const std::string& filename) {
        std::ifstream handle(filename.c_str());

        if(!handle.is_open()) {
            std::cerr << "Error: Could not find config file '" << filename << "'\n";
            std::exit(1);
        }

        std::string line;

        while(std::getline(handle, line)) {
            if (line.empty() || line[0] == '#') {
                continue;
            }

            std::size_t pos = line.find('=');

            if (pos == std::string::npos) // Switch
            {
                StringFormat::Strip(line);

                m_switches.emplace_back(line);
            }
            else
            {
                std::string key   = line.substr(0, pos++);
                std::string value = line.substr(pos, line.size() - pos);

                StringFormat::Strip(key);
                StringFormat::Strip(value);

                m_keys.emplace_back(key);
                m_values.emplace_back(value);
            }
        }

        handle.close();
    }

    void UpdateConfig(Config& config) {
        config.UpdateSwitches(m_switches);
        config.UpdateKeyValues(m_keys, m_values);
    }
};

#endif
