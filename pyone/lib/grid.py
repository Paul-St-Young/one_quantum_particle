import numpy as np

def read_grid_ang(filename):
    '''read data on regular grid in format:
    x y z data[x,y,z]

    return X,Y,Z,DATA
    '''

    # import data, get size
    grid_data = np.loadtxt(filename)
    x,y,z,dat = extract_grid_ang(grid_data)

    return x,y,z,dat

# end def read_grid_ang

def extract_grid_ang(grid_data):
    '''read data on regular grid in format:
    x y z data[x,y,z]

    return X,Y,Z,DATA
    '''

    # get size
    ndata = len(grid_data)
    size = int( round( ndata**(1./3) ) )

    # extract data
    dat = grid_data[:,3]

    # extract x,y,z !! in bohr
    b2a = 0.529
    grid_data = grid_data.reshape(size,size,size,5)
    x    = grid_data[:,0,0,0]/b2a
    y    = grid_data[0,:,0,1]/b2a
    z    = grid_data[0,0,:,2]/b2a

    return x,y,z,dat

# end def read_grid_ang

def init_lab3d(X,Y,Z):
    nx = len(X); ny = len(Y); nz = len(Z);
    dx = X[1] - X[0]; dy = Y[1] - Y[0]; dz = Z[1] - Z[0];
    # 7-point space-centered stencil 3D laplacian matrix
    lap3d = np.zeros([nx,ny,nz,nx,ny,nz]);
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                lap3d[i,j,k,i,j,k]       -= 2./dx**2.+2./dy**2+2./dz**2
                if (i-1>=0):
                    lap3d[i-1,j,k,i,j,k] += 1./dx**2.
                # end
                if (i+1<nx):
                    lap3d[i+1,j,k,i,j,k] += 1./dx**2.
                # end
                if (j-1>=0):
                    lap3d[i,j-1,k,i,j,k] += 1./dy**2.
                # end
                if (j+1<nx):
                    lap3d[i,j+1,k,i,j,k] += 1./dy**2.
                # end
                if (k-1>=0):
                    lap3d[i,j,k-1,i,j,k] += 1./dz**2.
                # end
                if (k+1<nx):
                    lap3d[i,j,k+1,i,j,k] += 1./dz**2.
                # end
            # end
        # end
    # end
    lap3d = lap3d.reshape([nx*ny*nz,nx*ny*nz])
    return lap3d
# end def init_lab3d

class Potential_Grid: # <- Scalar_Grid <- Rectangular_Domain

    def __init__(self):
        # units are not used yet, just for book-keeping for now
        domain_unit = "B"
        energy_unit = "H"
    # end def __init__

    def read_grid_ang(self,filename):
        '''read data on regular grid in format:
        x y z data[x,y,z], where x,y,z are in angstromm data is in Hatree
        !! assume data is in lexigraphical order.
        '''
        grid_data = np.loadtxt(filename)
        # get size
        ndata = len(grid_data)
        size = int( round( ndata**(1./3) ) )
        
        # extract data
        self.dat = grid_data[:,3]
        
        # extract x,y,z !! in bohr
        b2a = 0.529
        grid_data = grid_data.reshape(size,size,size,5)
        self.x    = grid_data[:,0,0,0]/b2a
        self.y    = grid_data[0,:,0,1]/b2a
        self.z    = grid_data[0,0,:,2]/b2a
        self.nx   = len(self.x)
        self.ny   = len(self.y)
        self.nz   = len(self.z)
        self.dx   = self.x[1] - self.x[0]
        self.dy   = self.y[1] - self.y[0]
        self.dz   = self.z[1] - self.z[0]
        self.measure = self.dx*self.dy*self.dz

        return self.x,self.y,self.z,self.dat
    # end def 

    def init_lap3d(self):
        X = self.x; Y=self.y; Z=self.z;
        nx = len(X); ny = len(Y); nz = len(Z);
        dx = X[1] - X[0]; dy = Y[1] - Y[0]; dz = Z[1] - Z[0];
        # 7-point space-centered stencil 3D laplacian matrix
        lap3d = np.zeros([nx,ny,nz,nx,ny,nz]);
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    lap3d[i,j,k,i,j,k]       -= 2./dx**2.+2./dy**2+2./dz**2
                    if (i-1>=0):
                        lap3d[i-1,j,k,i,j,k] += 1./dx**2.
                    # end
                    if (i+1<nx):
                        lap3d[i+1,j,k,i,j,k] += 1./dx**2.
                    # end
                    if (j-1>=0):
                        lap3d[i,j-1,k,i,j,k] += 1./dy**2.
                    # end
                    if (j+1<nx):
                        lap3d[i,j+1,k,i,j,k] += 1./dy**2.
                    # end
                    if (k-1>=0):
                        lap3d[i,j,k-1,i,j,k] += 1./dz**2.
                    # end
                    if (k+1<nx):
                        lap3d[i,j,k+1,i,j,k] += 1./dz**2.
                    # end
                # end
            # end
        # end
        lap3d = lap3d.reshape([nx*ny*nz,nx*ny*nz])
        self.lap3d = lap3d
        return lap3d
    # end def init_lab3d

    def grid_rep(self,func):
        return np.array([[[func(self.x[i],self.y[j],self.z[k]) for k in range(self.nz)] 
            for j in range(self.ny)] for i in range(self.nx)]).reshape(self.nx*self.ny*self.nz)
    # end def grid_rep

# end class Potential_Grid
