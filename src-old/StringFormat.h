#ifndef STRING_FORMAT__H__
#define STRING_FORMAT__H__

#include <string>
#include <vector>
#include <iostream>

class StringFormat {
public:
    static void LStrip(std::string& str) {
        auto it = str.begin();
        while (std::isspace(*it)) {
            it = str.erase(it);
        }
    }

    static void RStrip(std::string& str) {
        auto it = std::prev(str.end());
        while (std::isspace(*it)) {
            it = std::prev(str.erase(it));
        }
    }

    static void Strip(std::string& str) {
        LStrip(str);
        RStrip(str);
    }

    static void Split(std::vector<std::string>& list, const std::string& line, const std::string& delimiter) {
        std::size_t first = 0, second;
        do {
            second = line.find(delimiter, first);
            std::string element = line.substr(first, second - first);
            if (element.empty()) {
                std::cerr << "Error: List was badly formated: " << line << "\n";
                std::exit(1);    
            }
            Strip(element);
            list.emplace_back(element);
            first = second + delimiter.size();
        } while (second != std::string::npos);
    }
};

#endif