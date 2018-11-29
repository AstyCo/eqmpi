#include "globals.hpp"

void Asserter(const char *file, int line)
{
    std::cerr << "ASSERT at FILE:" << file << " LINE:"<< line << std::endl;
    exit(1);
}

MapGridDimensions getGridDimensions()
{
    MapGridDimensions mts;

    mts.insert(128, GridDimensions(4, 4, 8));
    mts.insert(256, GridDimensions(8, 4, 8));
    mts.insert(512, GridDimensions(8, 8, 8));

    return mts;
}

ConnectionDirection toCD(int i)
{
    return static_cast<ConnectionDirection>(i);
}
