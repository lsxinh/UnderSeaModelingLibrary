/**
 * @file min_grid.h
 * N-dimensional field for numerical minimization.
 */
#pragma once

#include <usml/types/seq_vector.h>

namespace usml {
namespace types {

/**
 * N-dimensional field for numerical minimization. The field has a seq_vector
 * axes in each dimension, and the field is defined at each point in this grid.
 * The gridded values are used to search for the axis offsets that would
 * result in a global minimum in the interpolated field.  This algorithm assumes
 * that global minimum for the interpolated field is bound by the points
 * around the minimum in the gridded field.
 */
class min_grid {
public:

    /**
     * Create N-dimensional field for numerical minimization.
     *
     * @param   num_dims    Number of dimension for this minimization.
     */
    min_grid( size_t num_dims )
            : _num_dims(num_dims) {
        _axis = new seq_vector*[_num_dims];
    }

    /**
     * Cleanup memory owned by this class.
     */
    virtual ~min_grid() {
        for (size_t dim = 0; dim < _num_dims; ++dim) delete _axis[dim];
        delete[] _axis;
    }

    /**
     * Retrieves axis associated with one dimension of the minimization.
     *
     * @param   dim     Dimension for this axis.
     * @return          Reference to this axis.
     */
    seq_vector* axis( size_t dim ) {
        return _axis[dim];
    }

    /**
     * Defines axis associated with one dimension of the minimization.
     *
     * @param   dim     Dimension for this axis.
     * @param   axis    Reference to this axis.
     */
    void axis( size_t dim, const seq_vector *axis ) {
        _axis[dim] = axis->clone();
    }

    /** Number of dimensions in this minimization. */
    const size_t num_dims() {
        return _num_dims;
    }

    /**
     * Used by sub-classes to define the function to be minimized.
     * Both arguments are defined as a double[] whose size is at least
     * as big as the number of dimensions in this minimization.
     *
     * @param   location    Location at which field value is desired.
     * @param   derivative  If this is not null, the first derivative
     *                      of the field at this point will also be computed.
     * @return              Value of the field at this point.
     */
    virtual double evaluate( double* location, double* derivative = NULL ) = 0;

    /**
     * Used by sub-classes to define the location at which to start the search.
     * This algorithm assumes that global minimum for the interpolated field
     * is bound by the points around this location.
     *
     * @return  Index numbers for the axes values that represent the
     *          minimum in the gridded field.  Returned as an array of size_t
     *          values whose size is equal to the number of dimensions
     *          in this minimization.
     */
    virtual const size_t* initial_best() = 0;

private:

    /** Number of dimensions in this minimization. */
    const size_t _num_dims;

    /** Axis associated with each dimension of the minimization. */
    seq_vector** _axis;

};

} /* namespace types */
} /* namespace usml */
