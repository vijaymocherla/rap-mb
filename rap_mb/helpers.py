# importing some paraphernelia 
import numpy as np
from scipy.constants import c, epsilon_0, hbar

import matplotlib as mpl
mpl.rc('text', usetex=True)
mpl.rc('font', family='sans-serif', serif='computer modern')
import matplotlib.pyplot as plt
# f     : focal length of the lens [cm]
# wl_nm : wavelength of the light [nm]
# r     : radius of the laser beam [cm]
# z     : distance from focal point [cm]
# vx    : velocity of molecular beam [m/s]
#
def zr(w0, wl_nm):
    """Calculates Rayleigh range in cm.
    """
    return np.pi*w0**2/(wl_nm*1e-7)

def w(z, wl_nm, r, f=0.0, focused=False):
    """Calculates transverse beam waist (in cm) at point along z the direction of propagation of the feild. 
    """
    w0 = r
    if focused:
        if abs(f)>1e-16:
            w0 = f*wl_nm*1e-7/(r*np.pi)
        else: raise Exception("!!!ERROR: Given focal length invalid.")
    waist = w0*np.sqrt(1+(z/zr(w0, wl_nm))**2)
    return waist  

def R(z, w0, wl_nm):
    """Caclulates the radius of curvature of the wavefronts in cm.
    """
    R =  z + zr(w0, wl_nm)**2/z
    return R 

def sweep_rate(z, vx, params):
    """Calculates the total sweep rate experienced by the molecule.
    """
    wl_nm, r, f = params
    omega_b = 2*np.pi*c/(wl_nm*1e-9)
    w0 = (wl_nm*1e-7*f/(np.pi*r))
    Rz = R(z, w0, wl_nm)*1e-2
    phi = vx**2 * omega_b/Rz/c
    return phi

def transit_time(z, vx, params):
    """Calculates the time taken by the molecule to cross the laser beam in seconds.
    """
    wl_nm, r, f = params
    beam_waist = w(z, wl_nm, r, f, focused=True)*1e-2
    return 2*beam_waist/vx

def delta(t, z, vel_xyz, params):
    """Calculates the time-dependent detunning, using the sweep rate.
    """
    wl_nm, r, f = params
    vx, vy, vz = vel_xyz
    omega_b = 2*np.pi*c/(wl_nm*1e-9)
    w0 = (wl_nm*1e-7*f/(np.pi*r))
    Rz = R(z, w0, wl_nm)*1e-2
    delta = vz*omega_b/c -  vx**2 *omega_b*t/(Rz*c)
    return delta

def Exy(t, z, P, vel_xyz, params):
    """Calculates the transverse Electric field of the laser beam.
    """
    wl_nm, r, f = params
    vx, vy, vz = vel_xyz
    # convert r,wx from [cm] to [m]
    wx = w(z, wl_nm, r, f, focused=True)*1e-2 
    wy = w(z, wl_nm, r, f=0.0, focused=False)*1e-2
    E0 = np.sqrt(4*P/(np.pi*epsilon_0*c*wx*wy))
    E = E0*np.exp(-(vx*t/wx)**2 -(vy*t/wy)**2)
    return E


# sampling velocities for Monte Carlo Simulations

def sample_velocities(v, delta_v, divergence):
    """Randomly samples (vx,vy,vz) for a beam travelling in x direction and with angular divergence along y and z.
    """
    # vx is longitudnal velocity of molecular beam 
    # (x-direction in lab frame)
    v_fwhm = v*delta_v/100
    v0 = np.random.normal(v, v_fwhm/2)  
    # velocity along the laser beam 
    # (z-direction in lab frame)
    alpha = np.random.normal(0, divergence)*np.pi/180.0
    beta = np.random.normal(0, divergence)*np.pi/180.0
    gamma = np.random.normal(0, divergence)*np.pi/180.0 
    vx = v0*np.cos(alpha)
    vy = v0*np.sin(beta)
    vz = v0*np.sin(gamma)
    return vx, vy, vz

def plot_property_data(data, x_key, param_key, title='', 
                    legend_title='', label_format='4.2f', legend_loc=(.80, .50),
                    colors={},
                    xlabel='', xlim=[], 
                    figsize=(5,4), save_fig=False, filename=''):
    fig, ax = plt.subplots(figsize=figsize)
    x = data[x_key]
    param_data = data[param_key]
    if colors=={}:
        tableau_palette = ['blue', 'orange','green','red','purple','brown','pink','gray','olive','cyan']
        colors = {ki: tableau_palette[i] for i, ki in enumerate(param_data.keys())}
    for ki in param_data:
        label = ('{k:'+label_format+'}').format(k=ki)
        ax.plot(x, param_data[ki], label=label, color=colors[ki])
    ax.legend(bbox_to_anchor=legend_loc, frameon=True, title=legend_title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r'$\textbf{N}_{ex}$')
    ax.set_xlim(np.min(x), np.max(x))
    if xlim != []:
        ax.set_xlim(xlim[0], xlim[1])
    ax.set_title(title)
    if save_fig:
        fig.savefig(filename, dpi=300)
    return fig, ax

# from scipy.integrate import solve_ivp

# def solve_bloch_eq(psi0, dipole, P, zi, vel_xyz, params, 
#                    method='DOP853', atol=1e-6, rtol=1e-3, nstep=500,
#                    dense_output=True):
#     """Solver for Optical-Bloch equations.
#     """
#     vx, vy, vz = vel_xyz
#     wl_nm, r, f = params
#     T = transit_time(zi, vx, params)
#     delta_zi = lambda t: delta(t, zi, vel_xyz, params) 
#     omega_zi = lambda t: dipole*Exy(t, zi, P, vel_xyz, params)/hbar
#     def func(t, y):
#         H = np.array([[0., -delta_zi(t), 0.], 
#                       [delta_zi(t), 0., omega_zi(t)],
#                       [0., -omega_zi(t), 0.]])
#         # print('%6.2E,%6.2E,%6.2E'%(t, delta_zi(t), omega_zi(t)))
#         return np.dot(H,y)
#     t_eval = None
#     if dense_output:
#         t_eval = np.linspace(-2*T, 2*T, nstep)
#     res = solve_ivp(func, t_span=(-2.*T, 2*T), t_eval=t_eval,
#                     y0=psi0, method=method, 
#                     dense_output=dense_output, rtol=rtol, atol=atol)
#     if not res.success:
#         print(res.message)
#         raise Exception("!!!ERROR: Couldn't propagate Optical-Bloch Equations")
#     return res.t, res.y


# def monte_carlo(run_params):
#     try:
#         wl_nm = run_params['wl_nm']
#         r = run_params['r']
#         f = run_params['f']
#         P = run_params['P']
#         # Molecular beam parameters
#         v = run_params['v']
#         delta_v = run_params['delta_v']
#         divergence = run_params['divergence']
#         dipole = run_params['dipole']
#         N_samples = run_params['N_samples']
#         Nt = run_params['N_t']
#         zi = run_params['zi']
#         U0 = run_params['U0']
#         V0 = run_params['V0']
#         W0 = run_params['W0']
#     except:
#         raise Exception("!!! ERROR: Please check input given in `run_params` of `type: dict()`")
#     psi0 = [U0, V0, W0]
#     params = wl_nm, r, f
#     vel_xyz = np.empty((N_samples, 3), dtype=np.float64)
#     w = np.empty(N_samples, dtype=np.float64)
#     for i in range(N_samples):
#         vel_xyz = sample_velocities(v, delta_v, divergence)
#         t, psi = solve_bloch_eq(psi0, dipole, P, zi, vel_xyz, params, method='DOP853', nstep=Nt)
#         w[i] = psi[2,-1]
#     w_avg = np.sum(w)/N_samples
#     n1_avg = (w_avg + 1.0)/2.0
#     return n1_avg




