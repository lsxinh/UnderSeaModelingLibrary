/**
 * @file min_golden.h
 * N-dimensional algorithm for Golden Sector Search.
 */
#pragma once

#include <usml/types/min_algorithm.h>

namespace usml {
namespace types {

/**
 * N-dimensional algorithm for the Golden Sector Search.  The Golden Sector
 * Search is a form of iterative bisection that uses uneven sections at each
 * step. It does not require the calculation of field derivatives.
 *
 * The algorithm assumes that the global minima is located in the interval
 * [a,b] around the minimum of the gridded field.  The best guess so far,
 * called "g", is initialized to the location of the minimum of the gridded
 * field.  At each step in the iteration, the algorithm creates a probe point,
 * called "p" in the larger of the two intervals [a,g] or [g,b].  The probe
 * point is placed such that it divides the larger interval into new sections
 * whose lengths are related by the Golden Ratio.  At each step, the interval
 * is reduced such that the lowest value so far is inside the new interval.
 */
class min_golden : public min_algorithm {
public:

    /**
     * Create N-dimensional field for numerical minimization.
     *
     * @param   num_dims    Number of dimension for this minimization.
     */
    min_golden( min_grid& field ) ;

    /**
     * Cleanup memory owned by this class.
     */
    virtual ~min_golden() ;

    /**
     * Search for the axis locations that would result in a global minimum
     * in the interpolated field.
     *
     * @return  Location at which field value is minimized.  Returned as an
     *          array of double values whose size is equal to the number of
     *          dimensions in this minimization.
     */
    virtual const double* minimize() ;

    /**
     * Maximum number of iterations to use in this search.
     * Defaults to 40.
     */
    size_t iteration_max() const {
        return _iteration_max;
    }

    /**
     * Maximum number of iterations to use in this search.
     * Defaults to 40.
     */
    void iteration_max( size_t iteration_max ) {
        _iteration_max = iteration_max;
    }

    /**
     * Iteration exits when range of function values is less than this.
     * Defaults to 1e-4.
     */
    double func_tolerance() const {
        return _func_tolerance;
    }

    /**
     * Iteration exits when range of function values is less than this.
     * Defaults to 1e-4.
     */
    void func_tolerance( double func_tolerance ) {
        _func_tolerance = func_tolerance;
    }

    /**
     * Iteration exits when relative interval_tolerance is less than this.
     * Defaults to 1e-4.
     */
    double interval_tolerance() const {
        return _interval_tolerance;
    }

    /**
     * Iteration exits when relative interval_tolerance is less than this.
     * Defaults to 1e-4.
     */
    void interval_tolerance( double interval_tolerance ) {
        _interval_tolerance = interval_tolerance;
    }

private:

    /**
     * Smaller side of a golden ratio interval.
     * Set to a constant ( 3 - sqrt(5) ) / 2  value.
     */
    const double _golden ;

    /** Maximum number of iterations to use in this search. */
    size_t _iteration_max ;

    /** Iteration exits when range of function values is less than this. */
    double _func_tolerance ;

    /** Iteration exits when relative interval_tolerance is less than this. */
    double _interval_tolerance ;

    /** Location on left side of [a,b] interval. */
    double* _pos_a ;

    /** Field value at _pos_a. */
    double _func_a ;

    /** Location on right side of [a,b] interval. */
    double* _pos_b ;

    /** Field value at _pos_b. */
    double _func_b ;

    /** Location of best "guess so far" in middle of interval. */
    double* _pos_g ;

    /** Field value at _pos_g. */
    double _func_g ;

    /** Probe point. */
    double* _pos_p ;

    /** Field value at _pos_p. */
    double _func_p ;

    /** Width of the [a,b] interval. */
    double* _width ;
};

} /* namespace types */
} /* namespace usml */
