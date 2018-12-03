#ifndef GLOBALS_HPP
#define GLOBALS_HPP

#include <mpi.h> // global include mpi just becouse we can

#include <map>

#include <iostream>
#include <sstream>

#include <cmath>

#define SH(x) sinh(x) // different floating point precision *hf, *h, *hl
#define CH(x) cosh(x) // different floating point precision *hf, *h, *hl

#define SSTR( x ) static_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

#define MY_ASSERT(x) if (!(x)) Asserter(__FILE__, __LINE__);
#define CHECK_INDEX(id, first, size) MY_ASSERT((id) >= (first) && (id) < (size));

typedef unsigned int uint;

#ifdef FLOAT_P
typedef float real;
#   ifndef MPI_TYPE_REAL
#       define MPI_TYPE_REAL MPI_FLOAT
#   endif
#else
typedef double real;
#   ifndef MPI_TYPE_REAL
#       define MPI_TYPE_REAL MPI_FLOAT
#   endif
#endif

enum MPI_SENDRECV_TAGS
{
    TAG_BOUNDARY_ELEMENTS
};


void Asserter(const char *file, int line);

inline real phi(real x, real y, real z)
{
    return sin(x) * SH(y - M_PI) * sin(z);
}

// Δ = div grad
// Δ phi = (d/dx^2 + d/dy^2 + d/dz^2) phi
inline real div_grad_phi(real x, real y, real z)
{
    return -phi(x, y, z); // sh(x)' = ch(x), ch(x)' = sh(x)
}

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
ConnectionDirection toCD(int i);

#endif // GLOBALS_HPP
