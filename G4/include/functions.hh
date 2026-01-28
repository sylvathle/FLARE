#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>

std::vector<std::string> split(const std::string&, char );

bool to_bool(const std::string&);

std::map<std::string, std::pair<bool, bool>> getBodies(const std::string& );

void printMap(const std::map<std::string, std::pair<bool, bool>>&);
