#ifndef COMMAND_LINE_PARSER__H__
#define COMMAND_LINE_PARSER__H__

#include <string>
#include <vector>
#include <iostream>

#include "Config.h"

class CommandLineParser {
private:
    std::vector<std::string> m_keys;
    std::vector<std::string> m_values;
    std::vector<std::string> m_switches;

public:
    CommandLineParser(const int argc, const char* argv[]) {
        ReadCommandLine(argc, argv);
    }

    void ReadCommandLine(const int argc, const char* argv[]) {
        for (int i = 1; i < argc; ++i) {
            std::string arg(argv[i]);
            
            if (arg.size() < 3 || !(arg[0] == '-' && arg[1] == '-')) {
                std::cerr << "Error: The format of the an argument was invalid (" << arg << ")\n";
                std::exit(1);
            }
            
            std::size_t pos = arg.find('=');
            
            if (pos == std::string::npos) // Switch
            {
                m_switches.emplace_back(arg.substr(2, arg.size() - 2));
            }
            else // Key - value pair
            {
                std::string key   = arg.substr(0, pos++);
                std::string value = arg.substr(pos, arg.size() - pos);

                if (key == "--" || value.empty()) {
                    std::cerr << "Error: Key or value was missing from argument (" << arg << ")\n";
                    std::exit(1);
                }

                m_keys.emplace_back(key.substr(2, key.size() - 2));
                m_values.emplace_back(value);
            }
        }
    }

    std::string ConfigFileName() const {
        for (std::size_t i = 0; i < m_keys.size(); ++i) {
            if (m_keys[i] == "config-file") {
                return m_values[i];
            }
        }
        return "ua.config";
    }

    void UpdateConfig(Config& config) const {
        config.UpdateSwitches(m_switches);
        config.UpdateKeyValues(m_keys, m_values);
    }
};

#endif
