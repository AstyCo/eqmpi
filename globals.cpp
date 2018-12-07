#include "globals.hpp"

#include "utils.hpp"

#include <iostream>
#include <sstream>

#include <stdlib.h> // exit

void Asserter(const char *file, int line)
{
    cnode.error(SSTR("ASSERT at FILE:" << file
                     << " LINE:"<< line << std::endl));
    exit(1);
}

MapGridDimensions getGridDimensions()
{
    MapGridDimensions mts;

    mts.insert(std::make_pair(1, GridDimensions(1, 1, 1)));
    mts.insert(std::make_pair(2, GridDimensions(2, 1, 1)));
    mts.insert(std::make_pair(4, GridDimensions(2, 2, 1)));
    mts.insert(std::make_pair(8, GridDimensions(2, 2, 2)));
    mts.insert(std::make_pair(16, GridDimensions(4, 2, 2)));

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

std::string CDtoString(ConnectionDirection cdir)
{
    switch (cdir) {
    case DIR_X: return std::string("DIR_X");
    case DIR_MINUS_X: return std::string("DIR_MINUS_X");
    case DIR_Y: return std::string("DIR_Y");
    case DIR_MINUS_Y: return std::string("DIR_MINUS_Y");
    case DIR_Z: return std::string("DIR_Z");
    case DIR_MINUS_Z: return std::string("DIR_MINUS_Z");
    case DIR_Y_PERIOD_FIRST: return std::string("DIR_Y_PERIOD_FIRST");
    case DIR_Y_PERIOD_LAST: return std::string("DIR_Y_PERIOD_LAST");
    case DIR_SIZE: return std::string("DIR_SIZE");
    }
}
