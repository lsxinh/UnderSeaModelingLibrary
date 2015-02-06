/**
 * @file min_golden.cc
 * N-dimensional algorithm for Golden Sector Search.
 */
#include <usml/ublas/ublas.h>
#include <usml/types/min_golden.h>

using namespace usml::types;

#define DEBUG_MIN_GOLDEN

/**
 * Create N-dimensional field for numerical minimization.
 */
min_golden::min_golden( min_grid& field ) :
    min_algorithm(field),
    _golden( (3.0-sqrt(5)) / 2.0 ),
    _iteration_max(40),
    _func_tolerance(1e-4),
    _interval_tolerance(1e-4)
{
    // create memory for N-dimensional locations

    _pos_a = new double[_num_dims];
    _pos_b = new double[_num_dims];
    _pos_g = new double[_num_dims];
    _pos_p = new double[_num_dims];
    _width = new double[_num_dims];

    // initialize locations in each dimension

    const size_t* index = _field.initial_best() ;
    for ( size_t dim=0 ; dim < _num_dims ; ++dim ) {
        _pos_a[dim] = _axis[dim]->operator[]( index[dim]-1 ) ;
        _pos_g[dim] = _axis[dim]->operator[]( index[dim] ) ;
        _pos_p[dim] = _pos_g[dim] ;
        _pos_b[dim] = _axis[dim]->operator[]( index[dim]+1 ) ;
        _width[dim] = _pos_b[dim] - _pos_a[dim] ;
    }

    // initialize field values at each location

    _func_a = _field.evaluate(_pos_a) ;
    _func_g = _field.evaluate(_pos_g) ;
    _func_p = _func_g ;
    _func_b = _field.evaluate(_pos_b) ;
}

/**
 * Cleanup memory owned by this class.
 */
min_golden::~min_golden() {
    delete[] _pos_a ;
    delete[] _pos_b ;
    delete[] _pos_g ;
    delete[] _pos_p ;
    delete[] _width ;
}

/**
 * Search for the axis locations that would result in a global minimum
 * in the interpolated field.
 */
const double* min_golden::minimize() {
    for (size_t iteration = 0; iteration < _iteration_max; ++iteration) {

        double norm  = 0 ;
        for (size_t dim = 0; dim < 3; ++dim) {

            // create a new probe point into larger section

            const double width_ag = _pos_g[dim] - _pos_a[dim];
            const double width_gb = _pos_b[dim] - _pos_g[dim];
            double step ;
            if (width_ag < width_gb) {          // probe the right side
                step = _golden * width_gb;
            } else {                            // probe the left side
                step = -_golden * width_ag;
            }
            _pos_p[dim] = _pos_g[dim] + step ;  // compute new probe point
            _func_p = _field.evaluate(_pos_p);  // eval function at probe point

            // display results for debugging

            #ifdef DEBUG_MIN_GOLDEN
                cout << "iteration=" << iteration << " dim=" << dim << endl
                     << "    _pos_a=" << _pos_a[dim] << " _func_a=" << _func_a
                     << "    _pos_g=" << _pos_g[dim] << " _func_g=" << _func_g << endl
                     << "    _pos_b=" << _pos_b[dim] << " _func_b=" << _func_b
                     << "    _pos_p=" << _pos_p[dim] << " _func_p=" << _func_p << endl ;
            #endif

            // decrease the size of the interval

            if (_func_p < _func_g) {                // >>> smallest at p <<<
                if (_pos_p[dim] < _pos_g[dim]) {    // new interval is a, p, g
                    _pos_b[dim] = _pos_g[dim];      // shorten interval on right
                    _func_b = _func_g;
                    _pos_g[dim] = _pos_p[dim];      // create new best guess
                    _func_g = _func_p;
                } else {                            // new interval is g, p, b
                    _pos_a[dim] = _pos_g[dim];      // shorten interval on left
                    _func_a = _func_g;
                    _pos_g[dim] = _pos_p[dim];      // create new best guess
                    _func_g = _func_p;
                }

            } else {                                // >>> smallest at g <<<
                if (_pos_p[dim] < _pos_g[dim]) {    // new interval is p, g, b
                    _pos_a[dim] = _pos_p[dim];      // shorten interval on left
                    _func_a = _func_p;
                } else {                            // new interval is a, g, p
                    _pos_b[dim] = _pos_p[dim];      // shorten interval on right
                    _func_b = _func_p;
                }
            }

            // compute the relative size of interval

            const double pos_diff = (_pos_b[dim] - _pos_a[dim]) / _width[dim];
            norm += pos_diff * pos_diff;

        } //loop over dim

        // exit if interval size is very small

        if ( norm < _interval_tolerance ) break ;

        // exit if function value has not changed very much

        const double _func_min = min(min(_func_a, _func_g), _func_b);
        const double _func_max = max(max(_func_a, _func_g), _func_b);
        if ( _func_max - _func_min < _func_tolerance ) break ;

    } // loop over iterations

    return _pos_g ;
}
