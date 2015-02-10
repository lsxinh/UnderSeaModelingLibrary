/**
 * @file offset_min_grid.cc
 * Compute wavefront offsets using numerical minimization
 */
#include <usml/waveq3d/offset_min_grid.h>

using namespace usml::waveq3d;

/**
 * Extract data from wavefront geometry and current target.
 */
offset_min_grid::offset_min_grid( size_t t1, size_t t2, size_t de, size_t az, const wave_queue& wave ) :
    min_grid(3), _target(*(wave.targets()), t1,t2)
{
    std::fill( _initial_best, _initial_best+3, 1 ) ;

    // construct 3-element axes for travel_time, source_de, source_az

    seq_linear travel_time(-wave.time_step(), wave.time_step(), 3);
    axis(0, &travel_time);

    double de_data[] = {
            wave.source_de(de - 1) - wave.source_de(de),
            0,
            wave.source_de(de + 1) - wave.source_de(de) };
    seq_data source_de(de_data, 3);
    axis(1, &source_de);

    double az_data[] = {
            wave.source_az(az - 1) - wave.source_az(az),
            0,
            wave.source_az(az + 1) - wave.source_az(az) };
    seq_data source_az(az_data, 3);
    axis(2, &source_az);

    // construct data grids for rho, theta, phi

    _wave_rho = new data_grid<double, 3>(axes());
    _wave_theta = new data_grid<double, 3>(axes());
    _wave_phi = new data_grid<double, 3>(axes());
    size_t index[3];
    for (size_t nd = 0; nd < 3; ++nd) {
        size_t d = de + nd - 1;
        index[1] = nd;
        for (size_t na = 0; na < 3; ++na) {
            size_t a = az + na - 1;
            index[2] = na;
            for (size_t nt = 0; nt < 3; ++nt) {
                wposition1* loc;
                switch (nt) {
                case 0:
                    loc = new wposition1(wave.prev()->position, d, a);
                    break;
                case 1:
                    loc = new wposition1(wave.curr()->position, d, a);
                    break;
                case 2:
                    loc = new wposition1(wave.next()->position, d, a);
                    break;
                }
                index[0] = nt;
                _wave_rho->data(index, loc->rho());
                _wave_theta->data(index, loc->theta());
                _wave_phi->data(index, loc->phi());
                delete loc;
            }
        }
    }
}

/**
 * Cleanup data created by this minimization.
 */
offset_min_grid::~offset_min_grid() {
    delete _wave_rho;
    delete _wave_theta;
    delete _wave_phi;
}

/**
 * Compute square of distance from point currently being evaluated
 * to the current acoustic target.
 */
double offset_min_grid::evaluate( double* location, double* derivative ) {
    if ( derivative ) {
        throw std::invalid_argument("derivative argument not supported");
    }
    _current.rho(_wave_rho->interpolate(location));
    _current.theta(_wave_theta->interpolate(location));
    _current.phi(_wave_phi->interpolate(location));
    return _current.distance2(_target);
}

/**
 * Compute distance from a location to the current acoustic target,
 * along each of the offset directions.
 */
double* offset_min_grid::distance( const double* location ) {
    double offset[3] ;
    std::fill( offset, offset+3, 0.0 ) ;
    for ( size_t n=0; n < 3 ; ++n ) {
        offset[n] = location[n] ;                   // offset along one axis
        _distance[n] = sqrt( evaluate(offset) );    // distance along that axis
        if (location[n] < 0.0) _distance[n]*=-1.0 ; // handle negative offsets
        offset[n] = 0.0 ;                           // reset offset
    }
    return _distance ;

}
