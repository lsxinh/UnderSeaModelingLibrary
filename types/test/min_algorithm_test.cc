/**
 * @example types/test/min_algorithm_test.cc
 */
#include <boost/test/unit_test.hpp>
#include <usml/types/types.h>
#include <stdexcept>

using namespace boost::unit_test;
using namespace usml::types;
using namespace usml::ublas;

/**
 * N-dimensional field for testing the ability to search for travel time,
 * launch DE, and launch AZ offsets given wavefront locations in the
 * neighborhood of the CPA. Test specifics are based on caustic path for
 * pedersen deep test with target at 3030 m. Scenario data was extracted
 * from wavefront *.nc file using Matlab, then hard coded into the
 * construction for this class.
 *
 * Builds a 3x3x3 grid of the rho, theta, and phi components of location
 * around the CPA.  Interpolates each of those grids to build locations
 * as a function travel time, launch DE, and launch AZ near CPA.  The
 * evaluate() method returns a measurement of the distance squared from each
 * of those locations to the current target location.
 */
class pedersen_caustic_scenario : public min_grid {

public:

    /**
     * Scenario data extracted from wavefront *.nc file using Matlab
     * for caustic path for pedersen deep test with target at 3030 m.
     */
    pedersen_caustic_scenario() :
        min_grid(3),
        _travel_time( 2.86, 0.01, 3 ),
        _source_de( 51.0, 0.2, 3 ),
        _source_az( -1.0, 1.0, 3 ),
        _target( 45.027219106612236, -45.0, -800.0 )
    {

        // wavefront locations in neighborhood of CPA

        double neighborhood[][3] = {
            {45.027225,-45.000672,-780.616843},
            {45.027230,-45.000000,-780.616843},
            {45.027225,-44.999328,-780.616843},
            {45.027209,-45.000672,-782.335240},
            {45.027213,-45.000000,-782.335240},
            {45.027209,-44.999328,-782.335240},
            {45.026005,-45.000642,-899.450746},
            {45.026009,-45.000000,-899.450746},
            {45.026005,-44.999358,-899.450746},
            {45.027289,-45.000674,-788.299896},
            {45.027293,-45.000000,-788.299896},
            {45.027289,-44.999326,-788.299896},
            {45.027271,-45.000674,-790.045127},
            {45.027275,-45.000000,-790.045127},
            {45.027271,-44.999326,-790.045127},
            {45.026063,-45.000644,-907.113985},
            {45.026067,-45.000000,-907.113985},
            {45.026063,-44.999356,-907.113985},
            {45.027351,-45.000675,-795.979509},
            {45.027356,-45.000000,-795.979509},
            {45.027351,-44.999325,-795.979509},
            {45.027334,-45.000675,-797.751246},
            {45.027338,-45.000000,-797.751246},
            {45.027334,-44.999325,-797.751246},
            {45.026121,-45.000645,-914.771588},
            {45.026125,-45.000000,-914.771588},
            {45.026121,-44.999355,-914.771588},
        };

        // construct a data grids from wavefront locations as a function of
        // travel_time, source_de, source_az

        std::fill( _initial_best, _initial_best+3, 0 ) ;
        axis( 0, &_travel_time ) ;
        axis( 1, &_source_de ) ;
        axis( 2, &_source_az ) ;
        _wave_rho = new data_grid<double, 3>(axes());
        _wave_theta = new data_grid<double, 3>(axes());
        _wave_phi = new data_grid<double, 3>(axes());
        size_t k = 0;
        size_t index[3];
        for (size_t nt = 0; nt < 3; ++nt) {
            index[0] = nt;
            for (size_t nd = 0; nd < 3; ++nd) {
                index[1] = nd;
                for (size_t na = 0; na < 3; ++na) {
                    index[2] = na;
                    double* data = neighborhood[k++];
                    wposition1 loc(data[0], data[1], data[2]);
                    _wave_rho->data(index, loc.rho());
                    _wave_theta->data(index, loc.theta());
                    _wave_phi->data(index, loc.phi());
                }
            }
        }
    }

    /**
     * Cleanup data created by this fixture.
     */
    ~pedersen_caustic_scenario() {
        delete _wave_rho;
        delete _wave_theta;
        delete _wave_phi;

    }

    /**
     * Compute square of distance from point currently being evaluated
     * to the current acoustic target.  This is the function to be minimized
     * by the search.
     *
     * @param location  Position in time, de angle, and az angle coordinates.
     * @return          Square of distance from this point to acoustic target.
     */
    virtual double evaluate( double* location, double* derivative = NULL ) {
        if ( derivative ) {
            throw std::invalid_argument("derivative argument not supported");
        }
        _current.rho(_wave_rho->interpolate(location));
        _current.theta(_wave_theta->interpolate(location));
        _current.phi(_wave_phi->interpolate(location));
        return _current.distance2(_target);
    }

    /**
     * Defines the location at which to start the search.
     *
     * @return  Index numbers (0,0,0) to indicate center of interval.
     */
    virtual const size_t* initial_best() {
        return _initial_best ;
    }

    /**
     * Travel times associated with wavefront locations
     * in neighborhood near CPA.
     */
    const seq_vector& travel_time() {
        return _travel_time;
    }

    /**
     * D/E angles associated with wavefront locations
     * in neighborhood near CPA.
     */
    const seq_vector& source_de() {
        return _source_de;
    }

    /**
     * AZ angles associated with wavefront locations
     * in neighborhood near CPA.
     */
    const seq_vector& source_az() {
        return _source_az;
    }

private:

    /**
     * Location at which to start the search.
     */
    size_t _initial_best[3] ;

    /**
     * Travel times associated with wavefront locations
     * in neighborhood near CPA.
     */
    seq_linear _travel_time;

    /**
     * D/E angles associated with wavefront locations
     * in neighborhood near CPA.
     */
    seq_linear _source_de;

    /**
     * AZ angles associated with wavefront locations
     * in neighborhood near CPA.
     */
    seq_linear _source_az;

    /**
     * Location of the point currently being evaluated.
     */
    wposition1 _current;

    /**
     * Location of the acoustic target.
     */
    wposition1 _target;

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

/**
 * Collection of tests in this suite
 */
BOOST_AUTO_TEST_SUITE(golden_search_test)

/**
 * @ingroup types_test
 * Test the ability of the golden_search algorithm to find wavefront offsets.
 */
BOOST_AUTO_TEST_CASE( golden_search ) {
    cout << "=== golden_search_test: golden_search ===" << endl;

    // construct N-dimensional algorithm for numerical minimization

    pedersen_caustic_scenario field ;
    min_golden search(field) ;
    search.iteration_max(2) ;

    // execute numerical minimization and display result

    const double* position = search.minimize() ;
    cout << "time=" << position[0] << "," << field.travel_time()[1]
         << "," << (position[0] - field.travel_time()[1]) << endl ;
    cout << "de=" << position[1] << "," << field.source_de()[1]
         << "," << (position[1] - field.source_de()[1]) << endl ;
    cout << "az=" << position[2] << "," << field.source_az()[1]
         << "," << (position[2] - field.source_az()[1]) << endl ;
}

/**
 * @ingroup types_test
 * Test the ability of the golden_search algorithm to find wavefront offsets.
 */
BOOST_AUTO_TEST_CASE( brent_search ) {
    cout << "=== golden_search_test: brent_search ===" << endl;
}

BOOST_AUTO_TEST_SUITE_END()
