import numpy as np
import scipy.optimize as op

# Zheng et al 2018 JCP 148:123329
# Function for solving distance from FRET efficiency using SAW-nu 

# input
# E: FRET efficiency; N: number of peptide bonds
# r0: R0 of FRET dye; prefactor: 0.55 nm.  We find 0.55 nm good for a well-mixed IDP sequence.
# r: distance vector for integration
# output
# d2saw: squared distance from SAW-nu
# nusaw: polymer scaling exponent from SAW-nu
def d2_saw(fe,N,r0,r,prefactor):
    rc=(r[1:]+r[:-1])/2.
    dr=r[1:]-r[:-1]
    res=fe.copy()
    new_nu=fe.copy()
    def func_fe(b,fe,N,prefactor):
        nu=np.log(np.sqrt(b)/prefactor)/np.log(N)
        prtmp=pr_saw(rc,b,nu)
        return np.sum(1./(1+(rc/r0)**6)*prtmp*dr)-fe
    for idx in range(len(fe)):
        ## (prefactor*N**0.5)**2 is the initial guess of squared distance
        res[idx]=op.fsolve(func_fe,x0=[(prefactor*N**0.5)**2],args=(fe[idx],N,prefactor))
    new_nu=np.log(np.sqrt(res)/prefactor)/np.log(N)
    return res,new_nu

def d2_saw_new(fe,N,r0,r,prefactor):
    rc=(r[1:]+r[:-1])/2.
    dr=r[1:]-r[:-1]
    res=fe.copy()
    new_nu=fe.copy()
    def func_fe(b,fe,N,prefactor):
        nu=np.log(np.sqrt(b)/prefactor)/np.log(N)
        prtmp=pr_saw(rc,b,nu)
        return (np.sum(1./(1+(rc/r0)**6)*prtmp*dr)-fe)**2
    for idx in range(len(fe)):
        ## (prefactor*N**0.5)**2 is the initial guess of squared distance
        res[idx]=op.minimize(func_fe,x0=[(prefactor*N**0.5)**2],args=(fe[idx],N,prefactor),method='Nelder-Mead').x
    new_nu=np.log(np.sqrt(res)/prefactor)/np.log(N)
    return res,new_nu

def d2_saw_new_r(fe,N,r0,r,prefactor):
    rc=(r[1:]+r[:-1])/2.
    dr=r[1:]-r[:-1]
    res=fe.copy()
    new_nu=fe.copy()
    def func_fe(b,fe,N,prefactor):
        nu=np.log(np.sqrt(b)/prefactor)/np.log(N)
        prtmp=pr_saw(rc,b,nu)
        return (1-np.sum(1./(1+(rc/r0)**6)*prtmp*dr)-fe)**2
    for idx in range(len(fe)):
        ## (prefactor*N**0.5)**2 is the initial guess of squared distance
        res[idx]=op.minimize(func_fe,x0=[(prefactor*N**0.5)**2],args=(fe[idx],N,prefactor),method='Nelder-Mead').x
    new_nu=np.log(np.sqrt(res)/prefactor)/np.log(N)
    return res,new_nu

def calc_pre(r):
    k=1.23e-32 * (1e8)**6   # angstrom^6/s^2
    tauc=5e-9  # s
    wh=0.85e9 # /s
    r2=10 # /s   #### now we have R2 from sPRE measurement!!!! It is sequence-dependent
    t=10.87e-3 #s
    r2sp=k/r**6 *(4*tauc+3*tauc/(1+wh**2*tauc**2))
    iox= r2*np.exp(-r2sp*t)/(r2+r2sp)
    return iox

def calc_pet(r):
    r0=1
    #beta=4
    return np.array(r<=r0,dtype=float)

def d2_saw_pet(pet,N,r,prefactor):
    rc=(r[1:]+r[:-1])/2.
    dr=r[1:]-r[:-1]
    res=pet.copy()
    new_nu=pet.copy()
    def func_pet(b,pet,N,prefactor):
        nu=np.log(np.sqrt(b)/prefactor)/np.log(N)
        prtmp=pr_saw(rc,b,nu)
        return (np.sum(calc_pet(rc)*prtmp*dr)-pet)**2
    for idx in range(len(pet)):
        ## (prefactor*N**0.5)**2 is the initial guess of squared distance
        res[idx]=op.minimize(func_pet,x0=[(prefactor*N**0.5)**2],args=(pet[idx],N,prefactor),method='Nelder-Mead').x
    new_nu=np.log(np.sqrt(res)/prefactor)/np.log(N)
    return res,new_nu

def d2_saw_pre(pre,N,r,prefactor):
    rc=(r[1:]+r[:-1])/2.
    dr=r[1:]-r[:-1]
    res=pre.copy()
    new_nu=pre.copy()
    def func_pre(b,pre,N,prefactor):
        nu=np.log(np.sqrt(b)/prefactor)/np.log(N)
        prtmp=pr_saw(rc,b,nu)
        return (np.sum(calc_pre(rc)*prtmp*dr)-pre)**2
    for idx in range(len(pre)):
        ## (prefactor*N**0.5)**2 is the initial guess of squared distance
        res[idx]=op.minimize(func_pre,x0=[(prefactor*N**0.5)**2],args=(pre[idx],N,prefactor),method='Nelder-Mead').x
    new_nu=np.log(np.sqrt(res)/prefactor)/np.log(N)
    return res,new_nu


weights=[0.0]*8
def d2_fpp(fe,pr,pe,N,r0,r,prefactor):
    rc=(r[1:]+r[:-1])/2.
    dr=r[1:]-r[:-1]
    res=pe.copy()
    new_nu=pe.copy()
    weight=fe.copy()
    def func_full(fits,fe,pr,pe,N,prefactor):
        b,w=fits
        nu=np.log(np.sqrt(b)/prefactor)/np.log(N)
        prtmp=pr_saw_gau(rc,w,b,nu,xc=0.8,sigma=0.45)
        return (np.sum(calc_pet(rc)*prtmp*dr)-pe)**2 + (np.sum(calc_pre(rc*10)*prtmp*dr)-pr)**2 + (np.sum(1./(1+(rc/r0)**6)*prtmp*dr)-fe)**2
    for idx in range(len(fe)):
        bounds=((0,(prefactor*N**0.6)**2),(0,1))
        minimum=op.minimize(func_full,x0=[(prefactor*N**0.5)**2,weights[idx]],args=(fe[idx],pr[idx],pe[idx],N,prefactor),bounds=bounds)
        res[idx],weight[idx]=minimum.x
        print(minimum.success)
        new_nu = np.log(np.sqrt(res)/prefactor)/np.log(N)
    return res, new_nu, weight


def d2_r0(fe,r01,pe,r02,N,r,prefactor):
    rc=(r[1:]+r[:-1])/2.
    dr=r[1:]-r[:-1]
    res=pe.copy()
    new_nu=pe.copy()
    weight=fe.copy()
    def func_full(fits,fe,r01,pe,f02,N,prefactor):
        b,w=fits
        nu=np.log(np.sqrt(b)/prefactor)/np.log(N)
        prtmp=pr_saw_gau(rc,w,b,nu,xc=1.3,sigma=0.41)
        return  (np.sum(1./(1+(rc/r01)**6)*prtmp*dr)-fe)**2 + (np.sum(1./(1+(rc/r02)**6)*prtmp*dr)-pe)**2
    for idx in range(len(fe)):
        bounds=((0,(prefactor*N**0.6)**2),(0,1))
        minimum=op.minimize(func_full,x0=[(prefactor*N**0.5)**2,weights[idx]],args=(fe[idx],r01,pe[idx],r02,N,prefactor),bounds=bounds)
        res[idx],weight[idx]=minimum.x
        print(minimum.success)
        new_nu = np.log(np.sqrt(res)/prefactor)/np.log(N)
    return res, new_nu, weight


def d2_rx(fe,r01,pe,r02,N,r,prefactor,x_c,sigma):
    rc=(r[1:]+r[:-1])/2.
    dr=r[1:]-r[:-1]
    res=pe.copy()
    new_nu=pe.copy()
    weight=fe.copy()
    def func_full(fits,fe,r01,pe,f02,N,prefactor):
        b,w=fits
        nu=np.log(np.sqrt(b)/prefactor)/np.log(N)
        prtmp=pr_saw_gau(rc,w,b,nu,xc=x_c,sigma=sigma)
        return  (np.sum(1./(1+(rc/r01)**6)*prtmp*dr)-fe)**2 + (np.sum(1./(1+(rc/r02)**6)*prtmp*dr)-pe)**2
    for idx in range(len(fe)):
        bounds=((0,(prefactor*N**0.6)**2),(0,1))
        minimum=op.minimize(func_full,x0=[(prefactor*N**0.5)**2,weights[idx]],args=(fe[idx],r01,pe[idx],r02,N,prefactor),bounds=bounds)
        res[idx],weight[idx]=minimum.x
        print(minimum.success)
        new_nu = np.log(np.sqrt(res)/prefactor)/np.log(N)
    return res, new_nu, weight


def d2_rsig(fe,r01,pe,r02,N,r,prefactor,sigma):
    rc=(r[1:]+r[:-1])/2.
    dr=r[1:]-r[:-1]
    res=pe.copy()
    new_nu=pe.copy()
    weight=fe.copy()
    def func_full(fits,fe,r01,pe,f02,N,prefactor):
        b,w=fits
        nu=np.log(np.sqrt(b)/prefactor)/np.log(N)
        prtmp=pr_saw_gau(rc,w,b,nu,xc=1.2,sigma=sigma)
        return  (np.sum(1./(1+(rc/r01)**6)*prtmp*dr)-fe)**2 + (np.sum(1./(1+(rc/r02)**6)*prtmp*dr)-pe)**2
    for idx in range(len(fe)):
        bounds=((0,(prefactor*N**0.6)**2),(0,1))
        minimum=op.minimize(func_full,x0=[(prefactor*N**0.5)**2,weights[idx]],args=(fe[idx],r01,pe[idx],r02,N,prefactor),bounds=bounds)
        res[idx],weight[idx]=minimum.x
        print(minimum.success)
        new_nu = np.log(np.sqrt(res)/prefactor)/np.log(N)
    return res, new_nu, weight


# P(r) of SAW-nu
def pr_saw(r,b,nu):
    g=(1.1615-1)/nu
    d=1./(1-nu)
    def tmp_saw(a,x):
        xc=(x[1:]+x[:-1])/2.
        dx=x[1:]-x[:-1]
        tmp=a[0]*xc**g*np.exp(-a[1]*xc**d)
        return (np.sum(tmp*4*np.pi*xc**2*dx)-1,np.sum(xc**2*tmp*4*np.pi*xc**2*dx)-1)
    x=r/np.sqrt(b)
    a=op.fsolve(tmp_saw,x0=[1.0,1.0],args=(x))
    prtmp=a[0]*x**g*np.exp(-a[1]*x**d)/(b)**1.5*4*np.pi*r**2
    return prtmp
# <R^2>/<Rg^2> as a function of the polymer scaling exponent
def ratio(nu):
    r=1.1615
    return 2*(r+2*nu)*(r+2*nu+1)/r/(r+1)

# Gaussian
# P(r)
def pr_gaussian(r,b):
    return 4*np.pi*r**2/(2./3*np.pi*b)**1.5*np.exp(-3*r**2/2/b)
def gaussian(x,xc,sigma):
    return 1/sigma/(2*np.pi)**0.5 * np.exp(-(x-xc)**2/2/sigma**2)
def pr_saw_gau(r,w,b,nu,xc,sigma):
        # two saw nu peaks at different nu with a weight between 0-1
        dr = r[1]-r[0]
        pr1 = (1-w)*pr_saw(r,b,nu)
        pr2 = w*gaussian(r,xc,sigma)
        return (pr1+pr2)*dr

# squared distance from FRET efficiency
def d2_gaussian(fe,r0,r):
    rc=(r[1:]+r[:-1])/2.
    dr=r[1:]-r[:-1]
    def func_fe(b,fe):
        prtmp=pr_gaussian(rc,b)
        return np.sum(1./(1+(rc/r0)**6)*prtmp*dr)-fe
    res=fe.copy()
    for idx in range(len(fe)):
        res[idx]=op.fsolve(func_fe,x0=[0.5],args=(fe[idx]))
    return res
