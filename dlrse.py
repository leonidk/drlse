from pylab import *
from skimage import filters, io
import scipy.ndimage as ndi
from time import sleep
#def filt_func(r, c, sigma = 1):
#   return np.exp(-np.hypot(r, c)/sigma)
#laplace = LPIFilter2D(filt_func)

def div(nx,ny):
    nxx,_ = np.gradient(nx)
    _,nyy = np.gradient(ny)
    return nxx + nyy

def delta(x, sigma):
    f = (0.5/sigma)*(1.0+cos(pi*x/sigma))
    b = (x<=sigma) & (x>=-sigma)
    f = f*b
    return f

def proc(img):
    #dx = filters.sobel_h(img)
    #dy = filters.sobel_v(img)
    #lap = filters.laplace(img)
    dx,dy = np.gradient(img)

def vNBounds(phi):
    phi[0,:] =  phi[1,:]
    phi[-1,:] = phi[-2,:]
    phi[:,0] =  phi[:,1]
    phi[:,-1] = phi[:,-2]

def distReg_p1(phi,curv):
    return 4*ndi.filters.laplace(phi)-curv
def distReg_p2(phi,dx,dy,mag):
    a = mag.copy()
    b = mag.copy()
    #a[a<0]=0
    #a[a>1]=0
    #a[a != 0] = 1

    #b[b <= 1] = 0
    #b[b > 1] = 1
    a = (mag >= 0) & (mag <= 1)
    b = (mag > 1)
    ps = a*sin(2.0*pi*mag)/(2*pi) + b*(mag-1)
    dps=((ps != 0)*ps + (ps == 0)) / ((mag != 0)*mag + (mag == 0))
    return div(dps*dx - dx, dps*dy -dy) + 4*ndi.filters.laplace(phi)
    #dps_den = mag.nonzero()
    #dps = ()/( (mag.nonzero())*()) #a*sin(2.0*pi*mag)/(2*pi) + b*(mag-1)

    #return 4*ndi.filters.laplace(phi)-curv

def drlse_edge(phi, edge, lambdap,mu,alpha,epsilon,timestep,iter_inner):
    vx, vy = np.gradient(edge)
    for i2 in range(iter_inner):
        vNBounds(phi) #edges are duplicated for no flux in or out of image
        dx,dy = np.gradient(phi)
        mag = np.sqrt((dx*dx)+(dy*dy))
        eps = 1e-3
        nx = dx/(mag+eps)
        ny = dy/(mag+eps)
        curv = div(nx,ny)
        #regTerm = distReg_p1(phi,curv)
        regTerm = distReg_p2(phi,dx,dy,mag)
        diracPhi = delta(phi,epsilon)
        areaTerm = diracPhi * mag
        edgeTerm = diracPhi * (vx*nx+vy*ny) + diracPhi*edge*curv
        phi += timestep*(mu*regTerm + lambdap*edgeTerm + alpha*areaTerm)
        #io.imshow(phi)
        #draw()
        #sleep(0.05)
        #show()
#params
def dslre(img):
    timestep = 1
    mu = 0.2/timestep
    iter_inner = 5
    iter_outer = 20
    iter_refine = 10
    lambdap = 5
    alpha = -3
    epsilon = 1.5
    sigma = 0.8

    smoothed = filters.gaussian_filter(img,sigma)
    dx,dy = np.gradient(smoothed)
    mag = (dx*dx)+(dy*dy)
    edge = 1.0/(1.0+mag)

    c0 = 2
    initialLSF = c0*np.ones(img.shape)
    initialLSF[25:35,20:25] -= c0
    initialLSF[25:35,40:50] -= c0
    phi = initialLSF
    for i in range(iter_outer):
        drlse_edge(phi,edge,lambdap,mu,alpha,epsilon,timestep,iter_inner)
    drlse_edge(phi,edge,lambdap,mu,0,epsilon,timestep,iter_refine)
    return phi

img = io.imread('gourd.bmp') #twocells
img = img.astype(np.float) 
ls = dslre(img)
imshow(ls)
show()
#io.imshow(dslre(img))
#show()
#io.imshow(img)
#show()