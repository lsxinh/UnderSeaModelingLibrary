/**
 * @file min_golden.h
 * N-dimensional algorithm for Golden Section Search.
 */
#pragma once

#include <usml/types/min_algorithm.h>

namespace usml {
namespace types {
/// @ingroup min_grid
/// @{

/**
 * N-dimensional algorithm for the Golden Section Search.  The Golden Section
 * Search is a form of iterative bisection that uses uneven sections at each
 * step. It does not require the calculation of field derivatives.
 *
 * The algorithm searches for a local minima in the interval [a,b]
 * around the minimum of the gridded field.  It assumes that only one
 * local minima exists in this interval. The best guess so far, called "g",
 * is initialized to the location of the minimum of the gridded
 * field.  At each step in the iteration, the algorithm creates a probe point,
 * called "p" in the larger of the two intervals [a,g] or [g,b].  The probe
 * point is placed such that it divides the larger interval into new sections
 * whose lengths are related by the Golden Ratio.  At each step, the interval
 * is reduced such that the lowest value so far is inside the new interval.
 *
 * In this implementation, several conditions can terminate the search:
 *
 *   - The iteration_max() property controls the maximum number of allowed
 *     iterations.  This defaults to 40.  Setting a hard limit prevents
 *     the algorithm from getting stuck in an infinite loop.
 *   - The iteration will also terminate if the difference between
 *     the minimum and maximum value of the field is less than the value
 *     defined by the func_tolerance() property.
 *   - The iteration will also terminate if the norm of the difference
 *     between the “a” and “b” limits of the interval, relative to the
 *     original size of the interval, is less than the value defined by
 *     the interval_tolerance() property.
 *
 * @xref Numerical Recipes: The Art of Scientific Computing,
 *       Third Edition (2007), Cambridge University Press, pp. 492-496.
 */
class min_golden : public min_algorithm {
public:

    /**
     * Initialize the a, b, g, and p locations, and their function values,
     * to setup the Golden Section Search.
     *
     * @param   field    N-dimensional field for numerical minimization.
     */
    min_golden( min_grid& field ) ;

    /**
     * Cleanup memory owned by this class.
     */
    virtual ~min_golden() ;

    /**
     * Iteratively places a the probe point such that it divides the larger
     * interval into new sections whose lengths are related by the
     * Golden Ratio.  At each iteration, the interval is reduced to
     * a smaller interval [a’,b’] such that the lowest value found
     * so far is inside the new interval.
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

/// @}
} /* namespace types */
} /* namespace usml */
