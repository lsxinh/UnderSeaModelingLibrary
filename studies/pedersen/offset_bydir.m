
clear ; close all ;

% wavefront step sizes in time, DE, AZ

delta = [ 0.01 0.2 1.0 ] ;

% distance squared values from caustic path
%   outer most loop is time
%   each group represents a DE angle
%   inner most loop is AZ angle

distance2 = zeros(3,3,3) ;
distance2(1,:,:) = [ [3174.05,377.063,3174.05]; [3107.78,312.547,3107.78]; [30690.5,28018.9,30690.5] ] ;
distance2(2,:,:) = [ [3007.52,204.063,3007.52]; [2940.17,138.484,2940.17]; [30587.9,27910.4,30587.9] ] ;
distance2(3,:,:) = [ [3056.77,246.844,3056.77]; [2988.16,180.063,2988.16]; [30692,28008.6,30692] ] ;

% compute gradient in all three directions, near center point

grad_t  = ( distance2(3,2,2) - distance2(1,2,2) ) / ( 2 * delta(1) ) ;
grad_de = ( distance2(2,3,2) - distance2(2,1,2) ) / ( 2 * delta(2) ) ;
grad_az = ( distance2(2,2,3) - distance2(2,2,1) ) / ( 2 * delta(3) ) ;
fprintf('grad: t=%f de=%f az=%f\n',grad_t,grad_de,grad_az);

[grad_de,grad_t,grad_az] = gradient( distance2, delta(2), delta(1), delta(3) ) ;
fprintf('grad: t=%f de=%f az=%f\n',grad_t(2,2,2),grad_de(2,2,2),grad_az(2,2,2));

% try alternate method for D/E gradient

left  = ( distance2(2,2,2) - distance2(1,2,2) ) / delta(1) ;
right = ( distance2(3,2,2) - distance2(2,2,2) ) / delta(1) ;
grad_t = ( left/sqrt(abs(left)) + right/sqrt(abs(right)) ) / ( 1/sqrt(abs(left)) + 1/sqrt(abs(right)) ) ;

left  = ( distance2(2,2,2) - distance2(2,1,2) ) / delta(2) ;
right = ( distance2(2,3,2) - distance2(2,2,2) ) / delta(2) ;
grad_de = distance2(2,2,2) * ( 1 / left + 1 / right ) / 2 ;
grad_de = 2 * ( 1 / left + 1 / right ) ;
grad_de = ( left + right ) / 2 ;
grad_de = ( left/sqrt(abs(left)) + right/sqrt(abs(right)) ) / ( 1/sqrt(abs(left)) + 1/sqrt(abs(right)) ) ;

left  = ( distance2(2,2,2) - distance2(2,2,1) ) / delta(3) ;
right = ( distance2(2,2,3) - distance2(2,2,2) ) / delta(3) ;
grad_az = distance2(2,2,2) * ( 1 / left + 1 / right ) / 2 ;
grad_az = ( left + right ) / 2 ;
grad_az = ( left/sqrt(abs(left)) + right/sqrt(abs(right)) ) / ( 1/sqrt(abs(left)) + 1/sqrt(abs(right)) ) ;

fprintf('grad: t=%f de=%f az=%f\n',grad_t,grad_de,grad_az);
