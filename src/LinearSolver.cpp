#include "Potentials.h"

#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

int main(int argc, const char* argv[]) {

	//std::cout << "Usage: <pmffile> <start>\n";

	std::ifstream handle(argv[1]);
	std::string line;
	std::vector<float> energy, distance;
	while (std::getline(handle, line)) {
		if (line.empty() || line[0] == '#') {
			continue;
		}
		std::size_t pos = line.find(',');
		if (pos == std::string::npos) // .dat file
		{
			distance.emplace_back(std::stod(line.substr(0, 8)));
			energy.emplace_back(std::stod(line.substr(8, 11)) * 0.40092); // kJ/mol -> kBT
		 }
		 else    // .csv file
		 {
			distance.emplace_back(std::stod(line.substr(0, pos++)));
			energy.emplace_back(std::stod(line.substr(pos, line.size() - pos)) * 0.40092);
		}
	}
	handle.close();

    Potential pot(energy, distance.front(), distance.back());

    const double startPosition = std::stod(argv[2]);

    double              zPosition       = startPosition;
    double              E;
    long double         area;
    long double         total_area      = 0.0;
    bool                isTouching      = false;
    const double        stepSize        = startPosition / 10000;

    while(true) {
	    E = pot.Value(zPosition, isTouching);
	    if(isTouching) {
		    break;
	    }
	    area = std::exp((long double) -1.0 * E);
	    // TODO: REMOVE ME ONCE DEBUGGING IS COMPLETE
	    if(std::isnan(area) || std::isinf(area)) {
		    std::cerr << "Warning: area was NaN\n";
		    break;
	    }
	    total_area += area;
	    zPosition -= stepSize;
    }

    std::cout << std::log(total_area * stepSize / startPosition) << "\n";

    return 0;
}
