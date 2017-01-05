#ifndef HELPER_H
#define HELPER_H

// #ifndef FEM_3D_BENCH
// #define IDEAL_IMPL
// // THIS IMPACTS THE PERFORMANCE SIGNIFICANTLY 
// #ifndef __INTEL_COMPILER
// #define COPY_SMART_EXPR
// #endif
// #endif

#include <Fastor.h>
#include <vector>
#include <fstream>
#include <sstream>

#ifndef _WIN32
#include <sys/resource.h>
#endif

using namespace Fastor;
using real = double;



#ifdef FEM_3D_BENCH

// compile time constants
#if POLYDEG == 1
constexpr int ndim          = 3;
constexpr int nnode         = 7234; 
constexpr int nelem         = 29998;
constexpr int nodeperelem   = 4;
constexpr int ngauss        = 4;
#elif POLYDEG == 2
constexpr int ndim          = 3;
constexpr int nnode         = 49264; 
constexpr int nelem         = 29998;
constexpr int nodeperelem   = 10;
constexpr int ngauss        = 5;
#elif POLYDEG == 3
constexpr int ndim          = 3;
constexpr int nnode         = 156026; 
constexpr int nelem         = 29998;
constexpr int nodeperelem   = 20;
constexpr int ngauss        = 14;
#elif POLYDEG == 4
constexpr int ndim          = 3;
constexpr int nnode         = 357602; 
constexpr int nelem         = 29998;
constexpr int nodeperelem   = 35;
constexpr int ngauss        = 31;
#elif POLYDEG == 5
constexpr int ndim          = 3;
constexpr int nnode         = 683906; 
constexpr int nelem         = 29998;
constexpr int nodeperelem   = 56;
constexpr int ngauss        = 53;
#elif POLYDEG == 6
constexpr int ndim          = 3;
constexpr int nnode         = 1165021; 
constexpr int nelem         = 29998;
constexpr int nodeperelem   = 84;
constexpr int ngauss        = 126;
#else
#error INTERPOLATION POLYNOMIAL DEGREE NOT UNDERSTOOD. SUPPLY -DPOLYDEG=<p> FLAG WHILE COMPILING
#endif

#else

// compile time constants
#if POLYDEG == 1
constexpr int ndim          = 2;
constexpr int nnode         = 165; 
constexpr int nelem         = 273;
constexpr int nodeperelem   = 3;
constexpr int ngauss        = 3;
#elif POLYDEG == 2
constexpr int ndim          = 2;
constexpr int nnode         = 603; 
constexpr int nelem         = 273;
constexpr int nodeperelem   = 6;
constexpr int ngauss        = 4;
#elif POLYDEG == 3
constexpr int ndim          = 2;
constexpr int nnode         = 1331; 
constexpr int nelem         = 273;
constexpr int nodeperelem   = 10;
constexpr int ngauss        = 7;
#elif POLYDEG == 4
constexpr int ndim          = 2;
constexpr int nnode         = 2298; 
constexpr int nelem         = 273;
constexpr int nodeperelem   = 15;
constexpr int ngauss        = 13;
#elif POLYDEG == 5
constexpr int ndim          = 2;
constexpr int nnode         = 3578; 
constexpr int nelem         = 273;
constexpr int nodeperelem   = 21;
constexpr int ngauss        = 19;
#elif POLYDEG == 6
constexpr int ndim          = 2;
constexpr int nnode         = 5107; 
constexpr int nelem         = 273;
constexpr int nodeperelem   = 28;
constexpr int ngauss        = 27;
#else
#error INTERPOLATION POLYNOMIAL DEGREE NOT UNDERSTOOD. SUPPLY -DPOLYDEG=<p> FLAG WHILE COMPILING
#endif

#endif


template<typename T>
T random() {
    return (T)rand()/RAND_MAX;
}



template<typename T, size_t rows, size_t cols>
Tensor<T,rows,cols> loadtxt(const std::string &filename)
{
    // Read to a Tensor

    T temp;
    std::ifstream datafile;
    datafile.open(filename.c_str());

    if(!datafile)
    {
        warn("Unable to read file");
    }

    Tensor<T,rows,cols> out_arr;

    for (int row=0; row<rows; ++row) {
        for (int col=0; col<cols; ++col) {
            datafile >> temp;
            out_arr(row,col) = temp;        
        }
    }

    datafile.close();

    return out_arr;
}



template<typename T, size_t M, size_t N, size_t O, size_t P, typename std::enable_if<M==N && M==2,bool>::type=0>
static inline Tensor<T,M,N,O,P> eye() {
    Tensor<T,M,N,O,P> out; out.zeros();
    out(0,0,0,0) = 1.0;
    out(0,0,1,1) = 1.0;
    out(1,1,0,0) = 1.0;
    out(1,1,1,1) = 1.0;
    return out;
}


template<typename T, size_t M, size_t N, size_t O, size_t P, typename std::enable_if<M==N && M==3,bool>::type=0>
static inline Tensor<T,M,N,O,P> eye() {
    Tensor<T,M,N,O,P> out; out.zeros();
    out(0,0,0,0) = 1.0;
    out(0,0,1,1) = 1.0;
    out(0,0,2,2) = 1.0;
    out(1,1,0,0) = 1.0;
    out(1,1,1,1) = 1.0;
    out(1,1,2,2) = 1.0;
    out(2,2,0,0) = 1.0;
    out(2,2,1,1) = 1.0;
    out(2,2,2,2) = 1.0;
    return out;
}


template<typename T, size_t M, size_t N, typename std::enable_if<M==N && M==3,bool>::type=0>
static inline Tensor<T,M,N> eye2() {
    Tensor<T,M,N> out; out.zeros();
    out(0,0) = 1.0;
    out(1,1) = 1.0;
    out(2,2) = 1.0;
    return out;
}

template<typename T, size_t M, size_t N, typename std::enable_if<M==N && M==2,bool>::type=0>
static inline Tensor<T,M,N> eye2() {
    Tensor<T,M,N> out; out.zeros();
    out(0,0) = 1.0;
    out(1,1) = 1.0;
    return out;
}


#endif // HELPER_H