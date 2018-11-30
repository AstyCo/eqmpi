#include "globals.hpp"

#include <iostream>
#include <sstream>

void Asserter(const char *file, int line)
{
    std::stringstream ss;
    ss << "ASSERT at FILE:" << file << " LINE:"<< line << std::endl;
    std::cerr << ss.str();
    exit(1);
}

MapGridDimensions getGridDimensions()
{
    MapGridDimensions mts;

    mts.insert(std::make_pair(32, GridDimensions(4, 4, 2)));
    mts.insert(std::make_pair(128, GridDimensions(4, 4, 8)));
    mts.insert(std::make_pair(256, GridDimensions(8, 4, 8)));
    mts.insert(std::make_pair(512, GridDimensions(8, 8, 8)));

    return mts;
}

ConnectionDirection toCD(int i)
{
    return static_cast<ConnectionDirection>(i);
}
