import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.time import Time
from astropy.coordinates import EarthLocation, SkyCoord
from PyAstronomy import pyasl
import tayph.system_parameters as sp
import astropy.units as u
import radvel
from PyAstronomy import modelSuite as ms
import seaborn as sns



plt.rc('font', size=12)          # controls default text sizes
plt.rc('axes', labelsize=18)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
plt.rc('ytick', labelsize=14)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
plt.rcParams["ytick.major.width"] = 2
plt.rcParams["xtick.major.width"] = 2
plt.rcParams["xtick.major.size"] = 5
plt.rcParams["ytick.major.size"] = 5
plt.rcParams["xtick.minor.size"] = 3.5
plt.rcParams["ytick.minor.size"] = 3.5
#print(plt.rcParams.keys())
plt.rcParams["ytick.minor.width"] = 1
plt.rcParams["xtick.minor.width"] = 1
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['legend.facecolor'] = 'white'
plt.rcParams['legend.edgecolor'] = 'white'
plt.rcParams['legend.framealpha'] = 1



class RVTraceEstimator:
    
    
    def __init__(self, dp, obs='paranal'):
        
        self.dp = dp
        self.obs = obs
        self.omega = None  # Define these variables in __init__ or set them elsewhere
        self.transitC = None
        self.period = None
        self.ecc = None
        self.bjds = None
        self.aRs = None
        self.Ms = None
        self.Mp = None
        self.transit = None
        self.RpRs = None
        self.T14 = None
        
        self.load_observation()
        self.read_orbital_configuration()
        self.calculate_planet_position()
        self.calculate_berv()
        self.RV_RM()
        self.RV_star()
        self.RV_planet()
        
    def load_observation(self):
        print('[INFO] Loading observation')
        
            
        # obstimes need to be in MJDs
        
        try:
            obstimes = ascii.read(f'{self.dp}/obs_times', comment="#")['col1']
        except:
            obstimes = []
        
        if len(obstimes) == 0: # so in case this file is empty, we need to produce them from scratch
            self.transitC = sp.paramget('Tc', self.dp)
            self.period = sp.paramget('P', self.dp)
            self.T14 = sp.paramget('T14', self.dp)
            
            obstimes = np.linspace(self.transitC-self.T14/2, self.transitC+self.T14/2, 1000)
            
        
        self.obstimes=obstimes
        mjds = Time(obstimes, format='mjd')
        self.bjds = mjds.jd
        
    
    def transit_ecc(self):
        print('[INFO] Determining in-transit expsoures')
        
        self.calculate_planet_position()
        
        rho   = np.sqrt(self.xp**2+self.yp**2) 
        dmin  = rho-self.RpRs #minimal distance between the planet ellipse and the origin
        dmax  = rho+self.RpRs #maximal distance between the planet ellipse and the origin
        
        
        transit = np.zeros_like(self.bjds)
        transit[self.zp < 0.] = 1. # planet is out of transit
        transit[dmin >= 1.] = 1 # planet is not overlapping with the stellar disk. 
        
        self.transit = transit        
    
    def read_orbital_configuration(self):
        print('[INFO] Reading orbital configuration')
        #self.load_observation()
        self.RA = sp.paramget('RA', self.dp)
        self.DEC = sp.paramget('DEC', self.dp)
        self.a = sp.paramget('a', self.dp)
        self.aRs = sp.paramget('aRstar', self.dp)
        self.orbinc = sp.paramget('inclination', self.dp)
        self.vsini = sp.paramget('vsini', self.dp)
        self.pob = sp.paramget('lampoo', self.dp)
        self.transitC = sp.paramget('Tc', self.dp)
        self.period = sp.paramget('P', self.dp)
        self.ecc = sp.paramget('ecc', self.dp)
        self.RpRs = sp.paramget('RpRstar', self.dp)
        self.K = sp.paramget('K', self.dp)
        self.vsys = sp.paramget('vsys', self.dp)
        
        
        if self.ecc != 0:
            
            self.omega = sp.paramget('omega', self.dp)  # Set the value of omega
            omega_bar = np.radians(self.omega)
            self.Ms = sp.paramget('Ms', self.dp) * u.Msun
            self.Mp = sp.paramget('Mp', self.dp) * u.Mjup
            self.T_per = radvel.orbit.timetrans_to_timeperi(self.transitC, self.period, self.ecc, omega_bar)
            
            
               
            
            self.transit_ecc()  
            #self.true_phase = 
            
            
        else:
            self.transit = sp.transit(self.dp)
            true_phases = sp.phase(self.dp)
            true_phases[true_phases > 0.5] -= 1. # making sure that they are indeed 
            
            self.true_phases = true_phases
    
        
    
    def calculate_planet_position(self):
        
        print('[INFO] Calculate planet position')
        
        if self.ecc==0:
            # circular orbit
            
            inclin_bar = np.radians(self.orbinc)
            
            self.true_anomaly = 2 * np.pi * self.true_phases
            x_pl = self.aRs * np.sin(self.true_anomaly)
            y_pl = -self.aRs * np.cos(self.true_anomaly) * np.cos(inclin_bar)
            z_pl = self.aRs * np.cos(self.true_anomaly) * np.sin(inclin_bar)
            node = np.radians(self.pob)
            
            
            self.xp = x_pl * np.cos(node) - y_pl * np.sin(node)
            self.yp = x_pl * np.sin(node) + y_pl * np.cos(node)
            self.zp = z_pl
            
        else:
            
            # degree to radians
            omega_bar = np.radians(self.omega)
            inclin_bar = np.radians(self.orbinc)
            

            ke = pyasl.KeplerEllipse(
                a=self.aRs,
                per=self.period,
                e=self.ecc,
                i=self.orbinc,  # in degrees
                Omega=0.,  # We don't care about the orientation in the sky
                w=self.omega,  # in degrees,
                tau=self.T_per
            )
            
            ecc_anomaly = ke.eccentricAnomaly(self.bjds)
            node = np.radians(self.pob)
            
            
            
            
            self.true_anomaly = 2. * np.arctan(np.sqrt((1. + self.ecc) / (1. - self.ecc)) * np.tan(ecc_anomaly / 2.))
            
            # in matrix form
            
            self.r0_p = self.aRs * np.array([
                                       (np.cos(ecc_anomaly) - self.ecc),
                                       np.sqrt(1. - self.ecc * self.ecc) * np.sin(ecc_anomaly),
                                        np.zeros_like(ecc_anomaly)
                                       ])
            
            
            rotation_angle = omega_bar #- np.pi/2
            r1_rotation_matrix =  np.array([
                                                  [np.sin(rotation_angle), np.cos(rotation_angle), 0],
                                                  [-np.cos(rotation_angle), np.sin(rotation_angle), 0],
                                                  [0,0,1]
                                       ])
                    
            
            # # create the eccentric orbit
            # self.X0_p = self.r0_p[0]
            # self.Y0_p = self.r0_p[1]
            
            self.r0_1 = r1_rotation_matrix @ self.r0_p
            
            self.X1_p = self.r0_1[0]
            self.Y1_p = self.r0_1[1]
            
            # # turn it w.r.t the argument of periastron
            # self.X1_p = self.X0_p * np.sin(omega_bar) + self.Y0_p * np.cos(omega_bar)
            # self.Y1_p = -self.X0_p * np.cos(omega_bar) + self.Y0_p * np.sin(omega_bar)
            
            # turn it w.r.t to the inclination
            
            r_pl_rotation_matrix = np.array([
                        [0, 1, 0],
                        [-np.cos(inclin_bar), 0, 0],
                        [np.sin(inclin_bar), 0, 0]
             
            ])
            
            self.r_pl = r_pl_rotation_matrix @ self.r0_1
            
            rp_rotation_matrix = np.array([
                        [np.cos(node), -np.sin(node), 0],
                        [np.sin(node), np.cos(node), 0],
                        [0, 0, 1]
            ])
            
            # x_pl = self.Y1_p
            # y_pl = -self.X1_p * np.cos(inclin_bar)
            # z_pl = self.X1_p * np.sin(inclin_bar)
            
            self.rp = rp_rotation_matrix @ self.r_pl
            
            self.xp = self.rp[0]
            self.yp = self.rp[1]
            self.zp = self.rp[2]
            
            # self.xp = x_pl * np.cos(node) - y_pl * np.sin(node)
            # self.yp = x_pl * np.sin(node) + y_pl * np.cos(node)
            # self.zp = z_pl
        

    def calculate_berv(self):
        print('[INFO] Calculating BERV')
        berv = []
        
    
        observatory = EarthLocation.of_site(self.obs)
        
        try:
            sc = SkyCoord(self.RA+' '+self.DEC, unit=(u.hourangle, u.deg))
        except:
            sc = SkyCoord(ra=self.RA * u.deg, dec=self.DEC*u.deg)
        for date in self.obstimes :
            barycorr = sc.radial_velocity_correction(obstime=Time(date,format='mjd'), location=observatory).to(u.km/u.s)
            berv.append(barycorr.value)
            
        self.berv = np.asarray(berv)
        
        
    def RV_RM(self):
        print('[INFO] Calculating the RV extent of the RM effect')
        self.calculate_planet_position()
        self.vstar = self.xp * self.vsini
        
        
    def RV_star(self):
        #calculate_planet_position(self)
        print('[INFO] Calculating the RV extent of the star')
        
        if self.ecc==0:
            self.RVstar = sp.RV_star(self.dp)
            
        else:
            rv = ms.KeplerRVModel()  # Changed from ms.KeplerRVModel()
            rv.assignValue({
                "per1": self.period,
                'K1': self.K,
                'e1': self.ecc,
                'tau1': self.T_per,
                'w1': self.omega,
                'mstar': self.Ms,
                'c0': 0.,
                'a1': self.a,
                'msini1': self.Mp * np.sin(np.radians(self.orbinc))
            })
            
            #print(self.bjds)
            
            RVs = rv.evaluate(self.bjds)  # evaluate the RV of the star at a given time
            self.RVstar = RVs
            
    def RV_planet(self):
        print('[INFO] Calculating the RV extent of the planet')
        #calculate_planet_position(self)
        
        if self.ecc==0:
            self.RVplanet = sp.RV(self.dp)
            
        else:
            self.RV_star()  # Call RV_star method of the class
            self.RVplanet = -1. * (self.RVstar) * (self.Ms / self.Mp).decompose()
              
    def plot_doppler_trace(self, RF='system', RV_ext=100, save=True):
        #print('[INFO] Calculating the RV extent of the star')
        
        # plotting in different restframes
        
        if self.ecc!=0:
            true_anomaly_sorted_idx = np.argsort(self.true_anomaly)
            true_anomaly_sorted = np.take_along_axis(self.true_anomaly, true_anomaly_sorted_idx, axis=0)
            RVstar_sorted = np.take_along_axis(self.RVstar, true_anomaly_sorted_idx, axis=0)
            vstar_sorted = np.take_along_axis(self.vstar, true_anomaly_sorted_idx, axis=0)  + RVstar_sorted
            vplanet_sorted = np.take_along_axis(self.RVplanet, true_anomaly_sorted_idx, axis=0)
            berv_sorted = np.take_along_axis(self.berv, true_anomaly_sorted_idx, axis=0)
            tel = - self.vsys + berv_sorted
            self.true_phases = self.true_anomaly / (2 * np.pi)
            transit_sorted =  np.take_along_axis(self.transit, true_anomaly_sorted_idx, axis=0)
            self.obtimes = np.take_along_axis(self.obstimes, true_anomaly_sorted_idx, axis=0)
            
        else:
            true_anomaly_sorted = self.true_anomaly
            RVstar_sorted = self.RVstar
            vstar_sorted = self.vstar.T +  RVstar_sorted
            vplanet_sorted = self.RVplanet
            berv_sorted = self.berv
            tel = - self.vsys + berv_sorted
            true_phases = self.true_phases * 2 * np.pi
            transit_sorted =  self.transit
        
        
        if RF=='star':
            
            self.RVs_RM = vstar_sorted - RVstar_sorted
            self.RVs_planet = vplanet_sorted - RVstar_sorted
            self.RVs_star = RVstar_sorted - RVstar_sorted
            self.RVs_tel = tel - RVstar_sorted
            
            print(f"[INFO] Plot in stellar restframe")
            
        
        elif RF=='planet':
            self.RVs_RM = vstar_sorted - vplanet_sorted
            self.RVs_planet = vplanet_sorted - vplanet_sorted
            self.RVs_star = RVstar_sorted - vplanet_sorted
            self.RVs_tel = tel - vplanet_sorted
            
            print(f"[INFO] Plot in planetary restframe")
        
        
        elif RF=='obs':
            self.RVs_RM = vstar_sorted - tel
            self.RVs_planet = vplanet_sorted - tel
            self.RVs_star = RVstar_sorted - tel
            self.RVs_tel = tel - tel
            
            print(f"[INFO] Plot in observatory restframe")
        
        elif RF=='berv':
      
            self.RVs_RM = vstar_sorted + self.vsys
            self.RVs_planet = vplanet_sorted + self.vsys
            self.RVs_star = RVstar_sorted + self.vsys
            self.RVs_tel = tel + self.vsys
            
            print(f"[INFO] Plot in berv restframe")

        elif RF=='system':
            self.RVs_RM = vstar_sorted
            self.RVs_planet = vplanet_sorted
            self.RVs_star = RVstar_sorted
            self.RVs_tel = tel

            print(f"[INFO] Plot in system restframe")
            
        
        
        elif 'all':
            
            print(f"[INFO] Plotting in all restframes")
            
            self.plot_doppler_trace(RF='star')
            self.plot_doppler_trace(RF='planet')
            self.plot_doppler_trace(RF='system')
            self.plot_doppler_trace(RF='obs')
            self.plot_doppler_trace(RF='berv')
            
        else:
            self.RVs_RM = vstar_sorted.T - RVstar_sorted
            self.RVs_planet = vplanet_sorted
            self.RVs_star = RVstar_sorted
            self.RVs_tel = tel_sorted

            print(f"[INFO] I don't know this restframe, I am plotting in the system restframe")


        
        #RV_ext = np.max(np.abs(np.concatenate()))
        
        if not RF=='all':
            plt.figure(figsize=(6,6))
            plt.title(f'Restframe {RF}')


            colors = sns.color_palette('colorblind', 4)

            # RM
            plt.plot(self.RVs_RM[transit_sorted != 1], true_anomaly_sorted[transit_sorted != 1] / (2 * np.pi), c=colors[0], label='RM', lw=3, ls='dashed')  # Use self.transit
            # planet
            plt.plot(self.RVs_planet, true_anomaly_sorted / (2 * np.pi),c=colors[1],  label='planet', lw=4)
            # star
            plt.plot(self.RVs_star, true_anomaly_sorted / (2 * np.pi), c=colors[2], ls='solid', label='star', lw=3)
            # tellurics:
            plt.plot(self.RVs_tel, true_anomaly_sorted / (2 * np.pi), lw=4, c=colors[3], ls='dotted', label='tellurics')  # Use self.vsys

            if transit_sorted[0] == 1:
                
                start_phase = true_anomaly_sorted[(transit_sorted != 1)][0] / (2 * np.pi)
                
                
                
                plt.axhline(start_phase, c='k', lw=2, ls='dashed')

            if transit_sorted[-1] == 1:
                end_phase = true_anomaly_sorted[(transit_sorted != 1)][-1] / (2 * np.pi)
                plt.axhline(end_phase, c='k', lw=2, ls='dashed')

            plt.xlim(-RV_ext, RV_ext)
            plt.ylabel('Phase')
            plt.xlabel('RV [km/s]')
            plt.legend(ncols=2)
            plt.tight_layout()
            

            #plt.show()
            if save:
                plt.savefig(f'{self.dp}RVTraceEstimator_{RF}.png')
