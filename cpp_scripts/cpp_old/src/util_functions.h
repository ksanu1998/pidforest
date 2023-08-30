#ifndef C_ISOLATION_UTIL_FUNCTIONS_H
#define C_ISOLATION_UTIL_FUNCTIONS_H

#include <iostream>
#include <cmath>
namespace utils
{
    inline void print__array(double * a, int size) {
        for (int i = 0; i < size; i++)
            std::cout<<a[i]<<" ";
        std::cout<<std::endl;
    }

    inline void print__array(const double * a, int size) {
        for (int i = 0; i < size; i++)
            std::cout<<a[i]<<" ";
        std::cout<<std::endl;
    }

    inline void print__array(int * a, int size) {
        for (int i = 0; i < size; i++)
            std::cout<<a[i]<<" ";
        std::cout<<std::endl;
    }

    inline double min(double a, double b) {
        if (a < b) return a;
        else return b;
    }

    inline bool isEqual(double x, double y)
    {
        const double epsilon = 1e-5; /* some small number such as 1e-5 */
        return std::abs(x - y) <= epsilon * std::abs(x);
        // see Knuth section 4.2.2 pages 217-218
    }
}
#endif
