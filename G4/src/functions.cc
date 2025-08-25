#include "functions.hh"


// Function to split a line by a delimiter
std::vector<std::string> split(const std::string& line, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(line);
    std::string token;

    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

// Function to convert a string to a boolean
bool to_bool(const std::string& str) {
    return str == "true" || str == "1";
}

// Function to load CSV and return a map
std::map<std::string, std::pair<bool, bool>> getBodies(const std::string& filename) {
    std::map<std::string, std::pair<bool, bool>> result;
    std::ifstream file(filename);
    std::string line;

    // Skip the header
    if (file.good()) {
        std::getline(file, line);
    }

    // Read the file line by line
    while (std::getline(file, line)) {
        std::vector<std::string> columns = split(line, ',');

        if (columns.size() == 3) {
            std::string body = columns[0];
            bool included = to_bool(columns[1]);
            bool detector = to_bool(columns[2]);

            result[body] = std::make_pair(included, detector);
        }
    }

    file.close();
    return result;
}

// Helper function to print the map for debugging
void printMap(const std::map<std::string, std::pair<bool, bool>>& data) {
    for (const auto& item : data) {
        std::cout << "Body: " << item.first
                  << " | Included: " << (item.second.first ? "true" : "false")
                  << " | Detector: " << (item.second.second ? "true" : "false")
                  << std::endl;
    }
}
