#ifndef NPTARGET_LIST__H__
#define NPTARGET_LIST__H__

#define BOOST_FILESYSTEM_VERSION_3
#define BOOST_FILESYSTEM_NO_DEPRECATED

#include <string>
#include <vector>
#include <iostream>
#include <boost/filesystem.hpp>

using Path = boost::filesystem::path;

class NPTargetList {
public:
    std::vector<std::string> m_paths;

public:
    NPTargetList(const std::vector<std::string>& targets) {
        if (targets.empty()) {
            std::cerr << "Error: No NP-target provided. Please provide a .np file or a directory containing .np files\n";
            std::exit(1);
        }
        
        for (const auto& target : targets) {
            CollectTargets(target);
        }
        
        if (m_paths.empty()) {
            std::cerr << "Error: Failed to find any NP files from provided target(s)\n";
            std::exit(1);
        }
    }

    void CollectTargets(const std::string& target) {
        Path targetPath(target);

        if (boost::filesystem::is_regular_file(targetPath) && targetPath.extension() == ".np") {
            m_paths.emplace_back(target);
        }
        else if (boost::filesystem::is_directory(targetPath)) {
            for (boost::filesystem::directory_iterator it{targetPath}; it != boost::filesystem::directory_iterator{}; ++it) {
                CollectTargets(it->path().string());
            }
        }
        else {
            std::cout << "Error: Unknown file type found! (" << target << ")\n";
        }
    }

    public:
        static std::string Filename(const std::string& target) {
            Path path(target);
            return path.filename().stem().string();
        }
};

#endif
