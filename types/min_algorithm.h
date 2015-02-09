/**
 * @file min_algorithm.h
 * N-dimensional algorithm for Golden Sector Search.
 */
#pragma once

#include <usml/types/min_grid.h>

namespace usml {
namespace types {
/// @ingroup min_grid
/// @{

/**
 * N-dimensional algorithm for numerical minimization.  Minimizing a field
 * requires the following programming tasks:
 *
 *   - Implement a subclass of min_grid that defines axes, implements
 *     the evaluate() method, and implements the initial_best() method.
 *   - Create an instance of this min_grid subclass.
 *   - Create an instance of one of the subclasses of min_algorithm.
 *   - Use the property setting accessors to adjust the algorithm's
 *     limits and tolerances (optional).
 *   - Use the minimize() method to Search for the axis locations that would
 *     result in a global minimum in the interpolated field.
 */
class min_algorithm {
public:

    /**
     * Create N-dimensional field for numerical minimization.
     *
     * @param   field    N-dimensional field for numerical minimization.
     */
    min_algorithm( min_grid& field ) : _field(field) {
    }

    /**
     * Cleanup memory owned by this class.
     */
    virtual ~min_algorithm() {
    }

    /**
     * Search for the axis locations that would result in a global minimum
     * in the interpolated field.
     *
     * @return  Location at which field value is minimized.  Returned as an
     *          array of double values whose size is equal to the number of
     *          dimensions in this minimization.
     */
    virtual const double* minimize() = 0 ;

protected:

    /** N-dimensional field for numerical minimization. */
    min_grid& _field;
};

/// @}
} /* namespace types */
} /* namespace usml */
