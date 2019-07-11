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


      //This is ugly but for testing purposes I want to see if it'll work before worrying about optimising the code.

//test to see if the amino acids are already defined in the config file
int useConfigAA = 0;

        for (std::size_t i = 0; i < m_keys.size(); ++i) {
if (m_keys[i] == "amino-acids") {
std::cout << "found amino acids in config file \n";
useConfigAA = 1;
}

}

//if not, we load these in from AAParams.dat
if(useConfigAA==0){
handle.open("AAParams.dat");

 if(!handle.is_open()) {
            std::cerr << "Error: Could not find AA parameter file \n";
            std::exit(1);
        }



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


    }

    void UpdateConfig(Config& config) {
        config.UpdateSwitches(m_switches);
        config.UpdateKeyValues(m_keys, m_values);
    }
};

#endif
