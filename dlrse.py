from pylab import *
from skimage import filters, io, color, morphology,exposure
from skimage.transform import hough_circle
from skimage.feature import peak_local_max, canny
from skimage.draw import circle_perimeter, circle
import scipy.ndimage as ndi
from time import sleep
import os, sys
from PIL import Image

def div(nx,ny):
    _,nxx = np.gradient(nx)
    nyy,_ = np.gradient(ny)
    return nxx + nyy

def delta(x, sigma):
    f = (0.5/sigma)*(1.0+np.cos(pi*x/sigma))
    b = (x<=sigma) & (x>=-sigma)
    return f*b

def vNBounds(phi):
    phi[0,:] =  phi[1,:]
    phi[-1,:] = phi[-2,:]
    phi[:,0] =  phi[:,1]
    phi[:,-1] = phi[:,-2]

def distReg_p1(phi,curv):
    return ndi.filters.laplace(phi)-curv

def distReg_p2(phi,dx,dy,mag):
    #dy, dx = np.gradient(phi)
    #mag = np.sqrt(dx**2+dy**2)
    a = (mag >= 0.) & (mag <= 1.)
    b = (mag > 1.)
    ps = a*np.sin(2.0*np.pi*mag)/(2.0*np.pi) + b*(mag-1.0)
    dps=((ps != 0.)*ps + (ps == 0.) ) / ((mag != 0.)*mag + (mag == 0.))
    return div(dps*dx - dx, dps*dy -dy) + ndi.filters.laplace(phi)


def drlse_edge(phi, edge, lambdap,mu,alpha,epsilon,timestep,iter_inner):
    vy, vx = np.gradient(edge)
    for i2 in range(iter_inner):
        vNBounds(phi) #edges are duplicated for no flux in or out of image
        dy,dx = np.gradient(phi)
        mag = np.sqrt((dx**2)+(dy**2))
        eps = 1e-6
        nx = dx/(mag+eps)
        ny = dy/(mag+eps)
        curv = div(nx,ny)
        
        #regTerm = distReg_p1(phi,curv)
        regTerm = distReg_p2(phi,dx,dy,mag)

        diracPhi = delta(phi,epsilon)
        #print nx.min(),nx.max(),curv.min(),curv.max(),regTerm.min(),regTerm.max(),diracPhi.min(),diracPhi.max()

        areaTerm = diracPhi * edge
        edgeTerm = diracPhi * (vx*nx+vy*ny) + diracPhi*edge*curv
        phi +=  timestep*(mu*regTerm + lambdap*edgeTerm + alpha*areaTerm)

#params
def dslre(img):
    timestep = 1.0
    mu = 0.2/timestep
    iter_basic = 1000
    iter_refine = 10
    lambdap = 5
    alpha = 1.5 # -3
    epsilon = 1.5
    sigma = 1.5

    smoothed = filters.gaussian_filter(img,sigma)
    dy,dx = np.gradient(smoothed)
    mag = (dx**2)+(dy**2)
    edge = 1.0/(1.0+mag)

    c0 = 2
    initialLSF = c0*np.ones(img.shape)
    initialLSF[10:50,10:50] = -c0

    #initialLSF[10:55,10:75] = -c0
    
    #initialLSF[25:35,20:25] -= c0
    #initialLSF[25:35,40:50] -= c0
    phi = initialLSF
    drlse_edge(phi,edge,lambdap,mu,alpha,epsilon,timestep,iter_basic)
    drlse_edge(phi,edge,lambdap,mu,0,epsilon,timestep,iter_refine)
    return phi
if False:
    img = io.imread('../Cell_03.png') #twocells.bmp
    if len(img.shape) == 3:
        img = img[:,:,1]
    img = img.astype(np.float) 
    phi = dslre(img)

    imshow(img)
    show()
    imshow(phi)
    show()
    seg, seg_n= ndi.label(phi < 0)
    print seg_n
    imshow(seg)
    show()
else:
    timestep = 1.0
    mu = 0.2/timestep
    iter_basic = 1000
    iter_refine = 0
    lambdap = 5
    alpha = 1.0#-3 #1.5 # -3
    epsilon = 1.5
    sigma = 0.8
    c0 = 2
    elem = morphology.disk(5)
    iml = None
    imgs = []
    for idx,fn in enumerate(sorted(os.listdir('.'))):
        img = io.imread(fn)
        if len(img.shape) == 3:
            img = img[:,:,1]
            img = img.astype(np.float) 
        imgs.append(img)
    meanImg = np.mean( np.array(imgs), axis=0 )
    stdImg = np.std( np.array(imgs),axis=0)

    for idx,fn in enumerate(sorted(os.listdir('.'))):
        img = io.imread(fn)
        orig_img = img.copy()
        if len(img.shape) == 3:
            img = img[:,:,1]
            img = img.astype(np.float) 
        
        #imr = (img-img.min())/(img.max()-img.min())
        #img = 255*exposure.equalize_adapthist(imr, clip_limit=0.01)

        if idx ==0:#True or idx == 0:
            initialLSF = c0*np.ones(img.shape)
            if True:
                initialLSF[10:50,10:50] = -c0 #cell 3
                #initialLSF[40:120,40:120] = -c0 #cell 8

                #initialLSF[10:-10,10:-10] = -c0
            else:
                edges = canny(img, sigma=3, low_threshold=10, high_threshold=50)
                #img = edges
                hough_radii = np.arange(20, 30, 3)
                hough_res = hough_circle(edges, hough_radii)

                centers = []
                accums = []
                radii = []

                for radius, h in zip(hough_radii, hough_res):
                    # For each radius, extract two circles
                    num_peaks = 1
                    peaks = peak_local_max(h, num_peaks=num_peaks)
                    centers.extend(peaks)
                    accums.extend(h[peaks[:, 0], peaks[:, 1]])
                    radii.extend([radius] * num_peaks)
                for idx in np.argsort(accums)[::-1]:
                    try:
                        center_x, center_y = centers[idx]
                        radius = radii[idx]
                        cx, cy = circle(center_y, center_x, radius)
                        initialLSF[cy, cx] = -c0
                    except:
                        pass
                initialLSF = morphology.erosion(initialLSF,elem)


            phi = initialLSF
        smoothed = filters.gaussian_filter(img,sigma)
        dy,dx = np.gradient(smoothed)
        mag = (dx**2)+(dy**2)
        edge = 1.0/(1.0+mag)

        drlse_edge(phi,edge,lambdap,mu,alpha,epsilon,timestep,iter_basic)
        drlse_edge(phi,edge,lambdap,mu,0,epsilon,timestep,iter_refine)
        
        if iml is None:
            iml = imshow(phi)
        else:
            iml.set_data(phi)
        #pause(0.01)
        draw()
        orig_img[:,:,0] = uint8(200*(phi < 0))
        #orig_img[:,:,1] += uint8(255*(phi < 0))
        #orig_img[:,:,2] += uint8(255*(phi < 0))

        pili = Image.fromarray(orig_img)
        pili.save('../out/' + fn)
        print fn
        #imshow(phi)
        #show()
        phi = morphology.erosion(phi,elem)
        phi[0,:] =  c0
        phi[-1,:] = c0
        phi[:,0] = c0
        phi[:,-1] = c0

        hough_trim = True
        if hough_trim:
            initialLSF = c0*np.ones(img.shape)

            edges = canny(img, sigma=3, low_threshold=10, high_threshold=50)
            #img = edges
            hough_radii = np.arange(20, 30, 3)
            hough_res = hough_circle(edges, hough_radii)

            centers = []
            accums = []
            radii = []

            for radius, h in zip(hough_radii, hough_res):
                # For each radius, extract two circles
                num_peaks = 1
                peaks = peak_local_max(h, num_peaks=num_peaks)
                centers.extend(peaks)
                accums.extend(h[peaks[:, 0], peaks[:, 1]])
                radii.extend([radius] * num_peaks)
            for idx in np.argsort(accums)[::-1]:
                try:
                    center_x, center_y = centers[idx]
                    radius = radii[idx]
                    cx, cy = circle(center_y, center_x, radius+3)
                    initialLSF[cy, cx] = -c0
                except:
                    pass
            phi = np.maximum(phi,initialLSF)
