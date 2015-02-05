%%
% brent_search.m - prototype use of Brent's Parabolic Interpolation 
% as a scheme for computing eigenray offsets.
%
% * Starts with a best guess g in the interval (a,b).
% * Subject to the condition that f(g) <= f(a) and f(g) <= f(b)
% * Creates a probe point p into largest interval (a,g) or (g,b).
% * Probe point is offset by golden fraction of interval width from g.
%
function brent_search
disp('*** brent_search ***');

% initialize scenario

threshold = 10^2 ;          % desired distance squared to target
golden = (3-sqrt(5))/2;     % golden offset from one edge of interval
earth_radius = 6371000.0 ;
d2r = pi/180 ;
target = define_target() ;
wavefront = load_neighborhood('pedersen_deep_raytrace_bydir.nc') ;
min_step = sqrt(eps) ;      % smallest realistic step size

% initialize offsets 

pos_a = [ -1 -1 -1 ] ;      % left side of interval
func_a = compute_distance2(pos_a) ;
pos_b = [ 1 1 1 ] ;         % right side of interval
func_b = compute_distance2(pos_b) ;
pos_g = [ 0 0 0 ] ;         % best guess so far
func_g = compute_distance2(pos_g) ;
pos_p = pos_g ;             % probe point
func_p = func_g ;

step=nan*ones(size(pos_p)); % step size for iteration
prev_step = step ;          % step size for previous iteration
past_step = step ;          % step size for two iterations ago
prev_fallback = false ;     % was fallback used last iteration?

% iterate solutions

for iteration=1:40
    if ( norm( pos_b - pos_a ) < 1e-4 || ...
            max([func_a func_b func_g]) - min([func_a func_b func_g]) < 1e-4 )
        break;
    end ;
    for dim = 1:3
        past_step(dim) = prev_step(dim) ;
        prev_step(dim) = step(dim) ;
        
        % compute commonly used differences
        % oriented such that widths are alway positive
        
        width_ab = pos_b(dim) - pos_a(dim) ;    % position diff a -> b
        width_ag = pos_g(dim) - pos_a(dim) ;    % position diff a -> g
        width_gb = pos_b(dim) - pos_g(dim) ;    % position diff g -> b
        func_ab = func_b - func_a ;             % func diff a -> b
        func_ag = func_g - func_a ;             % func diff a -> g
        func_gb = func_b - func_g ;             % func diff g -> b

        % estimate step size using inverse interpolation 
        
        r = width_gb * func_ag ;                % 1st term in demoninator
        s = width_ag * func_gb ;                % 2nd term in demoninator
        p = width_ag * r + width_gb * s ;       % parabolic numerator
        q = 2.0 * ( r - s ) ;                   % parabolic denominator
        
        if ( abs(q) > abs(p) )
            method = 'parabolic';               % inverse parabolic interp
            step(dim) = p / q ;
        elseif ( abs(func_ab) > abs(width_ab) )
            method = 'linear';                  % inverse linear interp
            step(dim) = - func_g * width_ab / func_ab ;
        else
            step(dim) = nan ;                   % fallback
        end
        
        % reject interpolation if
        %       1) step would jump outside of interval
        %       2) step size is smaller than min_step
        %       2) step size is larger than half old step
        % rejection is flagged by setting step size to NaN

        if ( ~isnan(step(dim)) )
            pos_p(dim) = pos_g(dim) + step(dim); % update probe location
            abs_step = abs(step(dim)) ;
            if ( pos_p(dim) <= pos_a(dim) || pos_p(dim) >= pos_b(dim) )
                method = 'outside of interval' ;
                step(dim) = nan ;       % jumped out of bounds
            elseif ( abs_step < min_step ) 
                method = 'step too small' ;
                step(dim) = nan ;       % step size too small
            else
                if ( prev_fallback )    % look at size of prev step
                    delta = abs( prev_step(dim) ) ;
                else                    % look at size of step before that
                    delta = abs( past_step(dim) ) ;
                end
                if ( delta > 0 && abs_step >= 0.5 * delta )
                    method = 'step too big' ;
                    step(dim) = nan ;   % step is too big
                end
            end
        end
        
        % fall back to golden section search if interpolation not ok
        
        prev_fallback = false ;
        if ( isnan(step(dim)) )
            if ( width_ag < width_gb )  % move probe the right side
                step(dim) = golden * width_gb ;
            else                        % move probe the left side
                step(dim) = - golden * width_ag ;
            end
            pos_p(dim) = pos_g(dim) + step(dim); % compute new probe point
            prev_fallback = true ;
        end
        
        func_p = compute_distance2(pos_p);  % eval function at probe point
        display();

        % decrease the size of the interval
        
        if ( func_p < func_g )              % >>> smallest at p <<<
            if ( pos_p(dim) < pos_g(dim) )  % new interval is a, p, g
                pos_b(dim) = pos_g(dim) ;   % shorten interval on right
                func_b = func_g ;
                pos_g(dim) = pos_p(dim) ;   % create new best guess
                func_g = func_p ;
            else                            % new interval is g, p, b
                pos_a(dim) = pos_g(dim) ;   % shorten interval on left
                func_a = func_g ;
                pos_g(dim) = pos_p(dim) ;   % create new best guess
                func_g = func_p ;
            end
        else                                % >>> smallest at g <<<
            if ( pos_p(dim) < pos_g(dim) )  % new interval is p, g, b
                pos_a(dim) = pos_p(dim) ;   % shorten interval on left
                func_a = func_p ;
            else                            % new interval is a, g, p
                pos_b(dim) = pos_p(dim) ;   % shorten interval on right
                func_b = func_p ;
            end
        end

    end     % loop over dim
end     % loop over iterations

return ;

    %%
    % Define scenario parameters for target location and CPA.
    % Based on caustic path for pedersen deep test with target at 3030 m.
    %
    function tgt = define_target()
        tgt.latitude = 45.027219106612236 ;
        tgt.longitude = -45.0 ;
        tgt.altitude = -800 ;
        tgt.rho = earth_radius + tgt.altitude ;
        tgt.theta = (90-tgt.latitude)*d2r ;
        tgt.phi = tgt.longitude*d2r ;
        tgt.de = 157 ;
        tgt.az = 5 ;
        tgt.time = 1+round(2.87/0.01) ;
    end

    %%
    % Load wavefront elements in neigborhood near CPA.
    %
    function wave = load_neighborhood(filename)
        wave = load_wavefront(filename) ;
        fprintf('time=%.3f source_de=%.3f source_az=%.3f\n',...
            wave.travel_time(target.time), ...
            wave.source_de(target.de), ...
            wave.source_az(target.az) ) ;
        offset = -1:1 ;
        t = target.time + offset ;
        d = target.de + offset ;
        a = target.az + offset ;
        
        wave.travel_time= wave.travel_time(t);
        wave.source_de  = wave.source_de(d);
        wave.source_az  = wave.source_az(a);
        wave.latitude   = wave.latitude(t,d,a);
        wave.longitude  = wave.longitude(t,d,a);
        wave.altitude   = wave.altitude(t,d,a);
        wave.surface    = wave.surface(t,d,a);
        wave.bottom     = wave.bottom(t,d,a);
        wave.caustic    = wave.caustic(t,d,a);
        wave.upper      = wave.upper(t,d,a);
        wave.lower      = wave.lower(t,d,a);
        wave.on_edge    = wave.on_edge(t,d,a);
        
        wave.rho        = earth_radius + wave.altitude ;
        wave.theta      = (90-wave.latitude)*d2r ;
        wave.phi        = wave.longitude*d2r ;
    end

    %%
    % Compute the distance squared for a particular combination of offsets
    % in the time, de, and az directions.  Each offset is assumed to be 
    % defined on the interval [-1,1]. Assumes linear interpolation on 
    % position from offsets. Relies on getting the "target" and 
    % "wavefront" variables from the parent function.  
    %
    function dist2 = compute_distance2( offset )
        interval = -1:1 ;
        point.rho   = interp3(interval,interval,interval,wavefront.rho,offset(1),offset(2),offset(3)) ;
        point.theta = interp3(interval,interval,interval,wavefront.theta,offset(1),offset(2),offset(3)) ;
        point.phi   = interp3(interval,interval,interval,wavefront.phi,offset(1),offset(2),offset(3)) ;
        dist2 = spherical_dist2( point, target ) ;
    end

    %%
    % Display a debugging print statement
    % 
    function display() 
        fprintf('iteration: %d dim: %d %s\n', iteration, dim, method ) ;
%         fprintf('\tpos_a=(%7.4f,%7.4f,%7.4f) func_a=%7.4f\n', ...
%             pos_a(1), pos_a(2), pos_a(3), func_a ) ;
%         fprintf('\tpos_g=(%7.4f,%7.4f,%7.4f) func_g=%7.4f\n', ...
%             pos_g(1), pos_g(2), pos_g(3), func_g ) ;
%         fprintf('\tpos_p=(%7.4f,%7.4f,%7.4f) func_p=%7.4f\n', ...
%             pos_p(1), pos_p(2), pos_p(3), func_p ) ;
%         fprintf('\tpos_b=(%7.4f,%7.4f,%7.4f) func_b=%7.4f\n', ...
%             pos_b(1), pos_b(2), pos_b(3), func_b ) ;
    end
end