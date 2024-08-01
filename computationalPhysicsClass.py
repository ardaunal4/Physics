import numpy as np
from scipy import integrate


class ComputationalPhysics:
    
    def D1Richardson(self, f, x, h):
        """
        This functional calculates first derivative 
        f: Function
        x: Argument of f
        h: Stepsize
        """
        return 1/(12*h) * ( f(x-2*h) - 8*f(x-h) + 8*f(x+h) - f(x+2*h) ) 
    
    def D2Richardson(self, f, x, h):
        """
        This functional calculates second derivative 
        f: Function
        x: Argument of f
        h: Stepsize
        """
        return 1/(12*h**2) * ( -f(x-2*h) + 16*f(x-h) - 30*f(x) + 16*f(x+h) - f(x+2*h) )

    def gradient(self, f, r, h):
        """
        This functional calculates divergences
        f: Function
        r: Position
        x: Argument of f
        h: Stepsize
        """
        x, y, z = r
        dx = self.D1Richardson(f, x, h)
        dy = self.D1Richardson(f, y, h)
        dz = self.D1Richardson(f, z, h)
        return np.array([dx, dy, dz])
    
    def divergence(self, f, r, h):
        """
        This functional calculates divergences
        f: Function
        r: Position
        x: Argument of f
        h: Stepsize
        """
        x, y, z = r
        dfdx = self.D1Richardson(f, x, h)
        dfdy = self.D1Richardson(f, y, h)
        dfdz = self.D1Richardson(f, z, h)
        return (dfdx + dfdy + dfdz)
    
    def curl(self, f, r, h):
        """
        This functional calculates divergences
        f: Function
        r: Position
        x: Argument of f
        h: Stepsize
        """
        x,y,z = r
        dfxdy = ( f(np.array([x, y+h, z]))[0] - f(np.array([x, y-h, z]))[0] ) / (2*h)
        dfxdz = ( f(np.array([x, y, z+h]))[0] - f(np.array([x, y, z-h]))[0] ) / (2*h)
        dfydx = ( f(np.array([x+h, y, z]))[1] - f(np.array([x-h, y, z]))[1] ) / (2*h)
        dfydz = ( f(np.array([x, y, z+h]))[1] - f(np.array([x, y, z-h]))[1] ) / (2*h)
        dfzdx = ( f(np.array([x+h, y, z]))[2] - f(np.array([x-h, y, z]))[2] ) / (2*h)
        dfzdy = ( f(np.array([x, y+h, z]))[2] - f(np.array([x, y-h, z]))[2] ) / (2*h)
        return np.array([ dfzdy-dfydz, dfxdz-dfzdx, dfydx-dfxdy ])
    
    def integralTrapezoidal(self, data):
        """
        This functional calculates integral by using Trapezoidal Method
        data: should be a numpy array with 2 arrays. 0th array consists of
        x values and 1st array consists of y values.
        """
        a = 0
        for i in range( len(data[0]) - 1 ):
            a = a + ( data[1, i+1] + data[1, i] ) / 2 * ( data[0, i+1] - data[0, i] )
        return a
    
    def integralTrapezoidalEQ(self, data):
        """
        This functional calculates integral by using equidistant Trapezoidal 
        Method
        data: should be a numpy array with 2 arrays. 0th array consists of
        x values and 1st array consists of y values.
        """
        return ( 1/2*data[1,0] + np.sum(data[1,1:-1]) + 1/2*data[1,-1] ) * (data[0,-1] - data[0,0]) / ( len(data[1]) - 1 )
    
    def integralSimpson(self, data):
        """
        This functional calculates integral by using Simpson Rule
        or Newton-Cortes equation.
        Note: There must be odd number of data points!
        data: should be a numpy array with 2 arrays. 0th array consists of
        x values and 1st array consists of y values.
        """
        return ( 1/3*data[1,0] + 4/3*np.sum(data[1,1:-1:2]) + 2/3*np.sum(data[1,2:-1:2]) + 1/3*data[1,-1] ) \
        * (data[0,-1] - data[0,0]) / ( len(data[1]) - 1 )
    
    def eulerODE(f, t0, y0, nmax, h):
        """
        f: Function with arguments arg1:t and arg2:y
        EXAMPLE FOR f: f = lambda t, y: y*t
        t0: Starting time
        y0: Starting value of y
        nmax: Number of iterations
        h: Stepsize
        """
        y = y0
        t = t0
        t_values = [t]
        y_values = [y]

        for i in range(1, nmax + 1):

            y = y + f(t, y) * h
            t = t + h
            t_values.append(t)
            y_values.append(y)

        return np.array([t_values, y_values])
    
    def eulerODE2(f, t0, y00, y10, nmax, h):
        """
        f: Function with arguments arg1:t and arg2:y arg2:y'
        EXAMPLE FOR f: f = lambda t, y, y1: y*t + y1 where
        y1 is the first derivative of the function
        t0: Starting time
        y00: Starting value of y(t)
        y10: Starting value of y'(t)
        nmax: Number of iterations
        h: Stepsize
        """
        y0 = y00
        y1 = y10
        t = t0
        t_values = [t]
        y0_values = [y0]
        y1_values = [y1]

        for i in range(1, nmax + 1):

            y0 = y0 + y1 * h
            y1 = y1 + f(t, y0, y1) * h
            t = t + h
            t_values.append(t)
            y0_values.append(y0)
            y1_values.append(y1)

        return np.array([t_values, y0_values, y1_values])
    
    def rk4(f, t0, y0, nmax, h):
        """
        Method: https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Classic_fourth-order_method
        f: Function
        t0: Starting time
        y0: Starting value of y
        nmax: Number of iterations
        h: Stepsize
        """
        y = y0
        t = t0
        t_values = [t]
        y_values = [y]

        for i in range(1, nmax + 1):

            k1 = h * f(t, y)
            k2 = h * f(t + h/2, y + k1/2)
            k3 = h * f(t + h/2, y + k2/2)
            k4 = h * f(t + h, y + k3)
            k = 1/6*k1 + 1/3*k2 +1/3*k3 + 1/6*k4
            y = y + k
            t = t + h
            t_values.append(t)
            y_values.append(y)
            
        return np.array([t_values, y_values])

    def rk45(f,t0,y0,nmax,h):
        """
        Method: https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Fehlberg
        f: Function
        t0: Starting time
        y0: Starting value of y
        nmax: Number of iterations
        h: Stepsize
        """
        y = y0
        t = t0
        t_values = [t]
        y_values = [y]

        for i in range(1, nmax + 1):

            k1 = h * f(t, y)
            k2 = h * f(t + h/4, y + k1/4)
            k3 = h * f(t + h*3/8, y + k1*3/32 + k2*9/32)
            k4 = h * f(t + h*12/13, y + k1*1932/2197 - k2*7200/2197 + k3*7296/2197)
            k5 = h * f(t + h, y + k1*439/216 - k2*8 + k3*3680/513 - k4*845/4104)
            k6 = h * f(t + h*1/2, y - k1*8/27 + k2*2 - k3*3544/2565 + k4*1859/4104 - k5*11/40)
            k = 16/135*k1 + 0*k2 +6656/12825*k3 + 28561/56430*k4 - 9/50*k5 + 2/55*k6
            y = y + k
            t = t + h
            t_values.append(t)
            y_values.append(y)
            
        return np.array([t_values, y_values])
    
    def scipyFirstOrderODESolver(f, ti, tf, starting_values, time_list, method = "RK45"):
        """
        f: Function with arguments arg1:t and arg2:y
        EXAMPLE FOR f: f = lambda t, y: y*t 
        ti: Initial time
        tf: Final time
        starting_values: Starting values of [y(t), y'(t), y"(t), ...]
        method: Solver method which could be
            1. RK45
            2. RK23
            3. DOP853
            4. Radau
            5. BDF
            6. LSODA
            Default is setted to RK45.
        time_list: Times at which to store the computed solution, 
        must be sorted and lie within ti-tf. If None (default), 
        use points selected by the solver. It must be a numpy array. 
        """
        return integrate.solve_ivp(f, [ti, tf], starting_values, method=method, t_eval = time_list)