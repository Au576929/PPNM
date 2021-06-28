This project is heavily based on the Runge-kutta12 implementaion made for the ODE homework, solving n-dimensional real-valued differential equations. 

The essential methodology of the ODE-solver remains the same, i.e. the differential-equation's rightside: f(x,y1,y2,...)
(note that solving n'th order differential equations can be done by e.q. y1=y(x), y2=dy/dx, y3=d²y/dx² ...) is evaluated, and the function is propagated through a step h in one and two steps,
the norm of difference between the two results is used to estimate the error of the given step. 

The first problem is a lot of variables need to have their type changed from e.g. double to double complex, using the complex.h library. Following that:

The two main issues of complex values for x and y are:
	1. x being complex means complex steps, which complicates things like testing whether "x+h>b", b being the end point of integration.
	
	2. The oútputs being complex means that the norms of the vectors needs to be calculated differently
	
the first problem is solved here by introducing a real variable t between 0 and 1, so that the stepper steps by dt*dh/dt, dh/dt being (b-a). So that when x starts in a when t goes to one x=b i.e. x(t)=a+(b-a)*t. 

The second problem is solved by introducing the norm function: |y|=sqrt(sum(|y_i|²)), which agrees with the traditional 2-norm on the real-axis. 

The driver takes as a variable a string, which is the name of a file to which the path is printed. Here the first coloumn is the realpart of x, second coloumn is the imaginary part of x,
the real part of y_1 followed by imaginary part of y_1 then similarly for y_2 and so on.


The plots are named after the harmonic-oscillator are described below:

SELF-EVALUATION:
As the updated algorithm solves ODE's of the required types, seemingly without issue, i would evaluate it at [10/10]. Thx ;-P

Harmonic_oscillator plots:
1 - solve from x=0 to x=20 with initial conditions y(0)=0 y'(0)=1 (y(x)=sin(x))
2 - solve from x=0 to x=5*I with initial conditions y(0)=1 y'(0)=0 (y(x)=cosh(imag(x)))
3 - solve from x=0 to x=5+5*I with inititial conditions y(0)=0 y'(0)=1 (y(x)=sin(x))


