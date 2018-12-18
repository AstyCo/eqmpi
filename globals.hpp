#ifndef GLOBALS_HPP
#define GLOBALS_HPP

#include <mpi.h> // global include mpi just becouse we can

#include <map>

#include <iostream>
#include <sstream>

#include <cmath>

#define SH(x) (sinh(x)) // different floating point precision *hf, *h, *hl
#define CH(x) (cosh(x)) // different floating point precision *hf, *h, *hl

#define SSTR( x ) static_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

#define MIN(a,b) ((a)<(b) ? (a) : (b))
#define MAX(a,b) ((a)<(b) ? (b) : (a))
#define ABS(a) ((a) > 0 ? (a) : (-(a)))

#define VAL_LY (double(M_PI))
#define VAL_LX (double(M_PI))
#define VAL_LZ (double(M_PI))
#define VAL_T (0.01)

#ifdef DEBUG
#define MY_ASSERT(x) if (!(x)) Asserter(__FILE__, __LINE__);
#define MY_ASSERT_X(x, text) if (!(x)) {cnode.error(text); Asserter(__FILE__, __LINE__);}
#define CHECK_INDEX(id, first, size) MY_ASSERT_X((id) >= (first) && (id) < (size), SSTR("Id " << id << " First " << first << " Size " << size));
#else
#define MY_ASSERT(x)
#define MY_ASSERT_X(x, text)
#define CHECK_INDEX(id, first, size)
#endif

typedef unsigned int uint;


#ifndef DOUBLE_P
typedef float real;
#   define MPI_TYPE_REAL MPI_FLOAT
#else
typedef double real;
#   define MPI_TYPE_REAL MPI_DOUBLE
#endif

void Asserter(const char *file, int line);

struct GridDimensions
{
    uint x;
    uint y;
    uint z;

    GridDimensions(uint x_ = 0, uint y_ = 0, uint z_ = 0)
        : x(x_), y(y_), z(z_)
    {}
};

typedef std::map<uint, GridDimensions> MapGridDimensions;

MapGridDimensions getGridDimensions();

enum ConnectionDirection
{
    DIR_X,
    DIR_MINUS_X,
    DIR_Y,
    DIR_MINUS_Y,
    DIR_Z,
    DIR_MINUS_Z,

    DIR_Y_PERIOD_FIRST,
    DIR_Y_PERIOD_LAST,

    DIR_SIZE
};

ConnectionDirection CDPair(ConnectionDirection cdir);
ConnectionDirection toCD(int i);
std::string CDtoString(ConnectionDirection cdir);

#endif // GLOBALS_HPP
