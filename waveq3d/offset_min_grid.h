/**
 * @file offset_min_grid.h
 * Compute wavefront offsets using numerical minimization
 */
#pragma once

#include <usml/waveq3d/wave_queue.h>

namespace usml {
namespace waveq3d {

using namespace usml::types;
using namespace usml::ublas;

/**
 * @internal
 * 3-dimensional field for computing wavefront offsets to a specific
 * target using numerical minimization.  Builds a 3x3x3 grid of the
 * rho, theta, and phi components of location around the CPA.
 * Interpolates each of those grids to build locations as a function
 * travel time, launch DE, and launch AZ near CPA.  The evaluate() method
 * returns a measurement of the distance squared from each
 * of those locations to the current target location.
 */
class offset_min_grid : public min_grid {

public:

    /**
     * Extract data from wavefront geometry and current target.
     */
    offset_min_grid( size_t t1, size_t t2, size_t de, size_t az, const wave_queue& wave ) ;

    /**
     * Cleanup data created by this minimization.
     */
    virtual ~offset_min_grid() ;

    /**
     * Defines the location at which to start the search.
     *
     * @return  Index numbers (0,0,0) to indicate center of interval.
     */
    virtual const size_t* initial_best() {
        return _initial_best ;
    }

    /**
     * Compute square of distance from point currently being evaluated
     * to the current acoustic target.  This is the function to be minimized
     * by the search.
     *
     * @param location  Position in time, de angle, and az angle coordinates.
     * @return          Square of distance from this point to acoustic target.
     */
    virtual double evaluate( double* location, double* derivative = NULL ) ;

    /**
     * Compute distance from a location to the current acoustic target,
     * along each of the offset directions.
     *
     * @param location  Position in time, de angle, and az angle coordinates.
     * @return          Square of distance from this point to acoustic target.
     */
    double* distance( const double* location ) ;

private:

    /**
     * Location of the acoustic target.
     */
    const wposition1 _target;

    /**
     * Location at which to start the search.
     */
    size_t _initial_best[3] ;

    /**
     * Location of the point currently being evaluated.
     */
    wposition1 _current;

    /**
     * Location of the point currently being evaluated.
     */
    double _distance[3] ;

    /**
     * Latitude component of wavefront location in the neighborhood of this target.
     */
    data_grid<double, 3>* _wave_rho;

    /**
     * Longitude component of wavefront location in the neighborhood of this target.
     */
    data_grid<double, 3>* _wave_theta;

    /**
     * Altitude component of wavefront location in the neighborhood of this target.
     */
    data_grid<double, 3>* _wave_phi;
};

}  // end of namespace waveq3d
}  // end of namespace usml
