/**
 * @example types/test/golden_search_test.cc
 */
#include <boost/test/unit_test.hpp>
#include <usml/types/types.h>

using namespace boost::unit_test;
using namespace usml::types;

/**
 * Test fixture to set up for a search in the neighborhood of
 * caustic path for pedersen deep test with target at 3030 m.
 * Scenario data was extracted from wavefront *.nc file using Matlab,
 * then hard coded into the construction for this class.
 */
class pedersen_caustic_scenario {

public:

    /**
     * Scenario data extracted from wavefront *.nc file using Matlab
     * for caustic path for pedersen deep test with target at 3030 m.
     */
    pedersen_caustic_scenario() :
        _travel_time( 2.86, 0.01, 3 ),
        _source_de( 51.0, 0.2, 3 ),
        _source_az( -1.0, 1.0, 3 ),
        _target( 45.027219106612236, -45.0, -800.0 )
    {
//        cout << "setup pedersen_caustic_scenario" << endl ;

        // Wavefront locations in neighborhood near CPA.

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

        seq_vector *axis[3] = { &_travel_time, &_source_de, &_source_az };
        _wave_rho = new data_grid<double, 3>(axis);
        _wave_theta = new data_grid<double, 3>(axis);
        _wave_phi = new data_grid<double, 3>(axis);
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
//        cout << "teardown pedersen_caustic_scenario" << endl;
        delete _wave_rho;
        delete _wave_theta;
        delete _wave_phi;

    }

    /**
     * Compute square of distance from point currently being evaluated
     * to acoustic target.  This is the function to be minimized
     * by the search.
     *
     * @param pos   Position in time, de angle, and az angle coordinates
     * @return      Square of distance from this point to acoustic target.
     */
    double compute_distance2( double pos[3] ) {
        _current.rho(_wave_rho->interpolate(pos));
        _current.theta(_wave_theta->interpolate(pos));
        _current.phi(_wave_phi->interpolate(pos));
        return _current.distance2(_target);
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
BOOST_FIXTURE_TEST_SUITE(golden_search_test, pedersen_caustic_scenario)

/**
 * @ingroup types_test
 * Test the ability of the golden_search algorithm to find wavefront offsets.
 */
BOOST_AUTO_TEST_CASE( golden_search ) {
    cout << "=== golden_search_test: golden_search ===" << endl;
    const double golden = (3.0-sqrt(5)) / 2.0  ;
    const size_t iteration_max = 40 ; // 40 ;
    const double func_tolerance = 1e-4 ;
    const double interval_tolerance = 1e-4 ;

    const double width[3] = {   // original width in each dimension
            travel_time()(2) - travel_time()(0),
            source_de()(2) - source_de()(0),
            source_az()(2) - source_az()(0),
    };

    // initialize point on left side of [a,b] interval

    double pos_a[3] = { travel_time()(0), source_de()(0), source_az()(0) } ;
    double func_a = compute_distance2( pos_a ) ;

    // initialize best "guess so far" as middle o interval

    double pos_g[3] = { travel_time()(1), source_de()(1), source_az()(1) } ;
    double func_g = compute_distance2( pos_g ) ;

    // initialize point on right side of [a,b] interval

    double pos_b[3] = { travel_time()(2), source_de()(2), source_az()(2) } ;
    double func_b = compute_distance2( pos_b ) ;

    // initialize "probe" point

    double pos_p[3] ;
    std::copy( pos_g, pos_g+3, pos_p ) ;
    double func_p = func_g ;

    // iterate solutions

    for (size_t iteration = 0; iteration < iteration_max; ++iteration) {

        double norm  = 0 ;
        for (size_t dim = 0; dim < 3; ++dim) {

            // create a new probe point into larger section

            const double width_ag = pos_g[dim] - pos_a[dim];
            const double width_gb = pos_b[dim] - pos_g[dim];
            double step ;
            if (width_ag < width_gb) {          // probe the right side
                step = golden * width_gb;
            } else {                            // probe the left side
                step = -golden * width_ag;
            }
            pos_p[dim] = pos_g[dim] + step ;    // compute new probe point
            func_p = compute_distance2(pos_p);  // eval function at probe point

            // display results for debugging

            cout << "iteration=" << iteration << " dim=" << dim << endl
                 << "    pos_a=" << pos_a[dim] << " func_a=" << func_a
                 << "    pos_g=" << pos_g[dim] << " func_g=" << func_g << endl
                 << "    pos_b=" << pos_b[dim] << " func_b=" << func_b
                 << "    pos_p=" << pos_p[dim] << " func_p=" << func_p << endl ;

            // decrease the size of the interval

            if (func_p < func_g) {              // >>> smallest at p <<<
                if (pos_p[dim] < pos_g[dim]) {  // new interval is a, p, g
                    pos_b[dim] = pos_g[dim];    // shorten interval on right
                    func_b = func_g;
                    pos_g[dim] = pos_p[dim];    // create new best guess
                    func_g = func_p;
                } else {                        // new interval is g, p, b
                    pos_a[dim] = pos_g[dim];    // shorten interval on left
                    func_a = func_g;
                    pos_g[dim] = pos_p[dim];    // create new best guess
                    func_g = func_p;
                }

            } else {                            // >>> smallest at g <<<
                if (pos_p[dim] < pos_g[dim]) {  // new interval is p, g, b
                    pos_a[dim] = pos_p[dim];    // shorten interval on left
                    func_a = func_p;
                } else {                        // new interval is a, g, p
                    pos_b[dim] = pos_p[dim];    // shorten interval on right
                    func_b = func_p;
                }
            }

            // compute the relative size of interval

            const double pos_diff = ( pos_b[dim] - pos_a[dim] ) / width[dim] ;
            norm += pos_diff * pos_diff ;

        } //loop over dim

        // exit if interval size is very small

        if ( norm < interval_tolerance ) break ;

        // exit if function value has not changed very much

        const double func_min = min(min(func_a, func_g), func_b);
        const double func_max = max(max(func_a, func_g), func_b);
        if ( func_max - func_min < func_tolerance ) break ;

    } // loop over iterations
}

/**
 * @ingroup types_test
 * Test the ability of the golden_search algorithm to find wavefront offsets.
 */
BOOST_AUTO_TEST_CASE( brent_search ) {
    cout << "=== golden_search_test: brent_search ===" << endl;
}

BOOST_AUTO_TEST_SUITE_END()
