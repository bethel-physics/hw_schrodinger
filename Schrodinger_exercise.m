%  schrodinger - Program to solve the Schrodinger equation 
%  for a free particle using the explicit FTCS, 
%  implicit FTCS, and Crank-Nicolson methods.
clear
close all

%% * Initialize parameters (grid spacing, time step, etc.)
i_imag = sqrt(-1);    % Imaginary i
N = 100;               % Number of grid points
L = 100;              % System extends from -L/2 to L/2
h = L/(N-1);          % Grid size
x = h*(0:N-1) - L/2;  % Coordinates  of grid points
h_bar = 1;  mass = 1; % Natural units
tau = input('Enter time step: ');
imethod=menu('Select method:','Explicit FTCS','Implicit FTCS','Crank-Nicolson');

%% * Set up the Hamiltonian operator matrix
ham = zeros(N);  % Set all elements to zero
coeff = _______;
for i=2:(N-1)
  ham(i,i-1) = _____;
  ham(i,i) = _____; %plus V(i) % Set interior rows
  ham(i,i+1) = _____;
end
% First and last rows for periodic boundary conditions
ham(1,N) = _____;   ham(1,1) = _____; ham(1,2) = _____;
ham(N,N-1) = _____; ham(N,N) = _____; ham(N,1) = _____;

if(imethod==1)     % Explicit FTCS matrix
    A = (eye(N) - i_imag*tau/h_bar*ham);
elseif(imethod==2) % Implicit FTCS
    A = _________________;
else               % Crank-Nicolson
    A = __________________;
end
			 
%% * Initialize the wavefunction 
x0 = 0;          % Location of the center of the wavepacket
velocity = 0.5;  % Average velocity of the packet
p0 = _____;
k0 = _____;       % Average wavenumber
sigma0 = L/10;   % Standard deviation of the wavefunction

Norm_psi = 1/(sqrt(sigma0*sqrt(pi)));  % Normalization
psi = Norm_psi * exp(i_imag*k0*x') .* ...
                      exp(-(x'-x0).^2/(2*sigma0^2));

%% * Plot the initial wavefunction
figure(1); clf;
plot(x,real(psi),'-',x,imag(psi),'--');
title('Initial wave function');
xlabel('x');  ylabel('\psi(x)'); legend('Real  ','Imag  ');
drawnow;  pause;

figure(2); clf;
theta=2*pi*(x+L/2)/L;
xt=cos(theta);
yt=sin(theta);
plot3(xt,yt,real(psi),'-',xt(1),yt(1),real(psi(1)),'ro');
hold on
plot3(xt,yt,imag(psi),'r-');
hold off
axis([-1 1 -1 1 -.3 0.3])
zlabel('Wave function');
view(-137.5,12);
drawnow; pause;

%% * Initialize loop and plot variables 
max_iter = L/(velocity*tau);      % Particle should circle system
plot_iter = max_iter/50;          % Produce 20 curves
p_plot(:,1) = ________;     % Record initial condition
iplot = 1;

figure(3); clf;
plot3(xt,yt,p_plot(:,1),'-',xt(1),yt(1),p_plot(1,1),'ro');
axis([-1 1 -1 1 -.02 0.07])
zlabel('Probability');
view(-137.5,12);
pause;

%% * Loop over desired number of steps (wave circles system once)
for iter=1:max_iter
	
  %* Compute new wave function based on the update matrix A
  psi = ______;  
  
  %* Periodically record values for plotting
  if( rem(iter,plot_iter) < 1 )   
    iplot = iplot+1;
    p_plot(:,iplot) = ________; 
    
    figure(2);
    plot3(xt,yt,real(psi),'-',xt(1),yt(1),real(psi(1)),'ro');
    hold on
    plot3(xt,yt,imag(psi),'r-');
    hold off
    axis([-1 1 -1 1 -.3 0.3])
    zlabel('Wave function');
    view(-137.5,12);
  
    figure(3);
    plot3(xt,yt,p_plot(:,iplot),'-',xt(1),yt(1),p_plot(1,iplot),'ro')
    axis([-1 1 -1 1 -.02 0.05])
    zlabel('Probability');
    view(-137.5,12);
    
    pause(.05)
  end

end

%% * Plot probability versus position at various times
figure(4); clf;
axisV = [-L/2 L/2 0 max(p_plot)]; % Fix axis min and max
pFinal = ________;
plot(x,p_plot(:,1:3:iplot),x,pFinal);    
xlabel('x'); ylabel('P(x,t)');
title('Probability density at various times');
   