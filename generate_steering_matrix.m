function A = generate_steering_matrix(M,d,DoA)

% d: Antenna spacing 

A = exp( -1i*2*pi*d*(0:M-1)'*sin(DoA*pi/180) );  % far-field model
