#include <string>
#include <fstream>
#include <iostream>

int main(const int argc, const char * argv[]) {
    
    if (argc != 2) {
        std::cout << "Usage: ./pdb2cg1 <file.pdb>\n";
        return 0;
    }
    
    std::ifstream handle(argv[1]);
    
    if (!handle.is_open()) {
        std::cerr << "Error: Could not find pdb file '" << argv[1] << "'\n";
        std::exit(1);
    }

    std::string line;

    while (std::getline(handle, line)) {
        if(line.size() > 3 && line.substr(0, 4) == "ATOM" && line.substr(13, 2) == "CA") {
            std::cout << line << "\n";
        }
        else if (line.size() > 5 && line.substr(0, 6) == "ENDMDL") {
            break;
        }
    }

    handle.close();
}
