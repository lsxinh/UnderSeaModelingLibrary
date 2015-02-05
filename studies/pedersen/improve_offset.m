%%
% improve_offset.m - prototype an improvement in compute_offset

clear ; close all

% conversion constants

earth_radius = 6371000.0 ;
d2r = pi/180 ;

% define target location and CPA parameters

target.latitude = 45.027219106612236 ;
target.longitude = -45.0 ;
target.altitude = -800 ;
target.rho = earth_radius + target.altitude ;
target.theta = (90-target.latitude)*d2r ;
target.phi = target.longitude*d2r ;
target.de = 157 ;
target.az = 5 ;
target.time = 1+round(2.87/0.01) ;

% load wavefront elements in this neigborhood

wavefront = load_netcdf_wavefront('pedersen_deep_raytrace_bydir.nc') ;
fprintf('time=%.3f source_de=%.3f source_az=%.3f\n',...
    wavefront.travel_time(target.time), ...
    wavefront.source_de(target.de), ...
    wavefront.source_az(target.az) ) ;
offset = -1:1 ;
t = target.time + offset ;
d = target.de + offset ;
a = target.az + offset ;
wavefront.travel_time = wavefront.travel_time(t);
wavefront.source_de = wavefront.source_de(d);
wavefront.source_az = wavefront.source_az(a);
wavefront.latitude = wavefront.latitude(a,d,t);
wavefront.longitude = wavefront.longitude(a,d,t);
wavefront.altitude = wavefront.altitude(a,d,t);
wavefront.surface = wavefront.surface(a,d,t);
wavefront.bottom = wavefront.bottom(a,d,t);
wavefront.caustic = wavefront.caustic(a,d,t);
wavefront.upper = wavefront.upper(a,d,t);
wavefront.lower = wavefront.lower(a,d,t);
wavefront.on_edge = wavefront.on_edge(a,d,t);

wavefront.rho = earth_radius + wavefront.altitude ;
wavefront.theta = (90-wavefront.latitude)*d2r ;
wavefront.phi = wavefront.longitude*d2r ;

n = (-1:1) ;
m = (-1:0.01:1);
points.rho = interp1( n, squeeze(wavefront.rho(2,2,:)), m ) ;
points.theta = interp1( n, squeeze(wavefront.theta(2,2,:)), m ) ;
points.phi = interp1( n, squeeze(wavefront.phi(2,2,:)), m ) ;

dist2 = spherical_dist2( points, target ) ;
dx = wavefront.travel_time(2) - wavefront.travel_time(1) ;
slope = (dist2(end)-dist2(1))/(2*dx)
figure ;
plot( m*dx, dist2) ;
grid
ylabel('Distance^2 (m^2)');
xlabel('Time Offset (sec)');

points.rho = interp1( n, squeeze(wavefront.rho(2,:,2)), m ) ;
points.theta = interp1( n, squeeze(wavefront.theta(2,:,2)), m ) ;
points.phi = interp1( n, squeeze(wavefront.phi(2,:,2)), m ) ;

dist2 = spherical_dist2( points, target ) ;
dx = wavefront.source_de(2) - wavefront.source_de(1) ;
slope = (dist2(end)-dist2(1))/(2*dx)
figure ;
plot( m*dx, dist2) ;
grid
ylabel('Distance^2 (m^2)');
xlabel('DE Offset (deg)');

points.rho = interp1( n, squeeze(wavefront.rho(:,2,2)), m ) ;
points.theta = interp1( n, squeeze(wavefront.theta(:,2,2)), m ) ;
points.phi = interp1( n, squeeze(wavefront.phi(:,2,2)), m ) ;

dist2 = spherical_dist2( points, target ) ;
dx = wavefront.source_az(2) - wavefront.source_az(1) ;
slope = (dist2(end)-dist2(1))/(2*dx)
figure ;
plot( m*dx, dist2) ;
grid
ylabel('Distance^2 (m^2)');
xlabel('AZ Offset (deg)');
