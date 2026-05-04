import numpy as np
import matplotlib.pyplot as plt
def potflow(params = None, xbounds = [-5, 5], ybounds= [-5, 5], levels=40, visualize = False):
    """
    Calculates the grid of psi, u, and v values of the flowfield 
    args:
        params: a 2D list. Each row of the list is considered to have the parameters of one elementary flow solution in the order [type, coordinate, strength]
        xbounds: the x bounds of the mesh grid on which to calculate the flowfield variables on 
        ybounds: the y bounds of the mesh grid on which to calculate the flowfield variables on 
        levels: the number of streamlines to plot
        visualize: to specify whether or not you want the function to automatically plot your flow visualization
    output: fig, ax, psival, U, V
        fig, ax: returns the figure containing the flow visualization to plot manually
        psival: 2D array of psi values across your grid 
        U: 2D array of x velocity values across your grid 
        V: 2D array of y velocity values across your grid 
    """
    if params is None:
        raise ValueError("params cannot be None")
    if len(params) == 0:
        raise ValueError("params cannot be empty")

    
    # Checking levels
    if not isinstance(levels, int) or levels <= 0:
        raise ValueError("levels must be a positive integer")

    
    def psi(X,Y,type, coord, strength):
        """
        Calculates the scalar field of the stream function of a single elementary flow solution over a grid of points
        Args:
            X,Y: 2D array of the x and y values of a grid of points
            type: a string that describes the type of elementary flow solution. It can be either "vor":vortex, "ufl": uniform flow, "ss": source/sink, "dbl": doublet
            coord: a 2 element array. For vor, ss, and dbl, coord is where it is centered and is written as [x0, y0]. for ufl, coord should be an array where the first element is the alpha in radians [alpha]
            strength: pass gamma for vortices, free stream velocity for ufl, lambda for source/sink, and kappa for doublet 
        """
        
        
        # Checking for each elementary flow type
        if not (isinstance(type, str)):
            raise TypeError("type argument must be a string")

        if type == "ufl":
            alpha = coord[0]
            psi_comp = strength*(Y*np.cos(alpha)-X*np.sin(alpha))
            return psi_comp
        else:
            x0 = coord[0]
            y0=coord[1]
            r = np.sqrt((X-x0)**2 + (Y-y0)**2)
            theta = np.atan2(Y-y0, X-x0)
        if type == "vor":
            psi_comp = strength*np.log(r)/(2*np.pi)
            return psi_comp
        elif type == "ss":
            psi_comp = strength*theta/(2*np.pi)
            return psi_comp
        elif type == "dbl":
            psi_comp = -strength*np.sin(theta)/(2*np.pi*r)
            return psi_comp
        else:
            raise ValueError("Invalid type argument")
    def vel(X,Y,type, coord, strength):
        """ 
        Returns an np array containing 2 2D arrays. The first 2D array is the u velocity of each point, the second 2D array is the v velocity at every point.
        takes the same arguments as the psi function
        """
        
        # Checks for each elementary flow solution type
        if type == "ufl":
            alpha = coord[0]
            u = strength*np.cos(alpha)
            v = strength*np.sin(alpha)
            U = u*np.ones_like(X)
            V = v*np.ones_like(X)
        # Other elementary flow solutions commonly use x0, y0 and r so they are precomputed here 
        else:
            x0 = coord[0]
            y0=coord[1]
            r = np.sqrt((X-x0)**2 + (Y-y0)**2)
        if type == "ss":
            U = (strength / (2 * np.pi)) * ((X-x0) / r**2)
            V = (strength / (2 * np.pi)) * ((Y-y0) / r**2)
        elif type == "vor":
            U = -(strength / (2 * np.pi)) * ((Y-y0) / r**2)
            V = (strength / (2 * np.pi)) * ((X-x0) / r**2)
        elif type == "dbl":
            U = -(strength / (2 * np.pi)) * ((X-x0)**2 - (Y-y0)**2) / r**4
            V = -(strength / (2 * np.pi)) * (2 * (X-x0) * (Y-y0)) / r**4
        return np.array([U,V])


    x = np.linspace(xbounds[0], xbounds[1], levels*10) #resolution of the grid is determined by the number of levels specified
    y = np.linspace(xbounds[0], ybounds[1], levels*10)
    X,Y = np.meshgrid(x,y)
    psival = np.zeros_like(X)
    uv = np.zeros_like(X)
    points = []


    #adding up the psi and velocity components of each elementray flow solution
    for param in params:
        if len(param) != 3:
            raise ValueError("Each param must have 3 elements: [type, coord, strength]")
        typeval, coordval, strengthval = param
        psival = psival + psi(X,Y,str(typeval), coordval, strengthval)
        uv = uv+vel(X,Y,str(typeval), coordval, strengthval)
        if typeval != "ufl": 
            print('noted')
            points.append([coordval[0], coordval[1]])

    points = np.array(points)
    vmag = np.linalg.norm(uv, axis=0)
    print(np.max(vmag))


    #creates the figure containing the flow visualization
    plt.ion()
    fig, ax = plt.subplots()
  
    ax.contour(X,Y, psival, levels=[0, 1], colors = "blue", linewidths=1.5, linestyles = 'solid')
    cmesh=ax.pcolormesh(X,Y, vmag, vmax=np.max(vmag), vmin=np.min(vmag), cmap='viridis', shading='gouraud')
    if len(points) != 0: ax.scatter(points[:,0], points[:,1], c='black')
    cbar = plt.colorbar(cmesh)
    cbar.set_label('Velocity magnitude (m/s)')

    #optionally plots the figure
    if visualize:
        fig.canvas.draw()
        plt.pause(0.1) 

        input("Press Enter to close...")
        plt.show()
    

    return fig, ax, psival, uv[0], uv[1]

