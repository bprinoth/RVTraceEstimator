import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.time import Time
from astropy.coordinates import EarthLocation, SkyCoord
from PyAstronomy import pyasl
#import tayph.system_parameters as sp
import astropy.units as u
import radvel
from PyAstronomy import modelSuite as ms



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
    
    
    def __init__(self, input_params, obsdate, Tdur, obs='paranal'):
        
        self.transitC, self.period, self.ecc, self.omega, self.a, self.aRs, self.orbinc, self.vsini, self.pob, self.Ms, self.Mp, self.RpRs, self.K, self.vsys, self.RA, self.DEC, self.n_exp  = input_params
        
        self.obs = obs
        self.T14 = float(Tdur)
        obsdate_utc = Time(obsdate,  format='isot', scale='utc') #+ 2400000.5
        self.obsdate = (obsdate_utc.tdb).jd

        # overwrite transitC with the value that comes after the transit start. So the first transit c after. 
        
        deltat = self.obsdate - self.transitC # time has passed since the transit c 
        #print(deltat / self.period)
        self.transitC = self.transitC + np.int64(deltat / self.period) * self.period

        self.load_observation()
        self.read_orbital_configuration()
        self.calculate_planet_position()
        self.calculate_berv()
        self.RV_RM()
        self.RV_star()
        self.RV_planet()

        
    def load_observation(self):
        print('[INFO] Generate observation')
        
        obstimes = np.linspace(self.obsdate, self.obsdate+(self.T14)/24, int(self.n_exp)) - 2400000.5
            
        #print(obstimes, self.transitC)
        self.obstimes=obstimes
        mjds = Time(self.obstimes, format='mjd')
        self.bjds = mjds.jd 
 
    def transit_ecc(self):
        print('[INFO] Determining in-transit expsoures')
        
        self.calculate_planet_position()
        
        rho   = np.sqrt(self.xp**2+self.yp**2) 
        #print(self.xp, self.yp, self.zp)
        dmin  = rho-self.RpRs #minimal distance between the planet ellipse and the origin
        dmax  = rho+self.RpRs #maximal distance between the planet ellipse and the origin
        
        transit = np.zeros_like(self.bjds)
        transit[self.zp < 0.] = 1. # planet is out of transit
        transit[dmin >= 1.] = 1 # planet is not overlapping with the stellar disk. 
        
        self.transit = transit     
        print(self.transit)   
    
    def read_orbital_configuration(self):
        print('[INFO] Reading orbital configuration')
        #self.load_observation()
        # self.RA = sp.paramget('RA', self.dp)
        # self.DEC = sp.paramget('DEC', self.dp)
        # self.a = sp.paramget('a', self.dp)
        # self.aRs = sp.paramget('aRstar', self.dp)
        # self.orbinc = sp.paramget('inclination', self.dp)
        # self.vsini = sp.paramget('vsini', self.dp)
        # self.pob = sp.paramget('lampoo', self.dp)
        # self.transitC = sp.paramget('Tc', self.dp)
        # self.period = sp.paramget('P', self.dp)
        # self.ecc = sp.paramget('ecc', self.dp)
        # self.RpRs = sp.paramget('RpRstar', self.dp)
        # self.K = sp.paramget('K', self.dp)
        # self.vsys = sp.paramget('vsys', self.dp)
        
        
        omega_bar = np.radians(self.omega)
        self.T_per = radvel.orbit.timetrans_to_timeperi(self.transitC, self.period, self.ecc, omega_bar)
        
        self.transit_ecc()  
    
    def calculate_planet_position(self):
            
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
        
        #print(ecc_anomaly, self.T_per)
        
        
        self.true_anomaly = 2. * np.arctan(np.sqrt((1. + self.ecc) / (1. - self.ecc)) * np.tan(ecc_anomaly / 2.))
        
        # in matrix form
        
        self.r0_p = self.aRs * np.array([
                                   (np.cos(ecc_anomaly) - self.ecc),
                                   np.sqrt(1. - self.ecc * self.ecc) * np.sin(ecc_anomaly),
                                    np.zeros_like(ecc_anomaly)
                                   ])
        
        
        rotation_angle = omega_bar - np.pi/2
        r1_rotation_matrix =  np.array([
                                              [np.cos(rotation_angle), -np.sin(rotation_angle), 0],
                                              [np.sin(rotation_angle), np.cos(rotation_angle), 0],
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

        
        self.RV_star()  # Call RV_star method of the class
        self.RVplanet = -1. * (self.RVstar) * (self.Ms * u.Msun / (self.Mp * u.Mjup)).decompose()
              
    def calculate_RV_traces(self, RF='system'):
        
        
        true_anomaly_sorted_idx = np.argsort(self.true_anomaly)
        self.true_anomaly_sorted = np.take_along_axis(self.true_anomaly, true_anomaly_sorted_idx, axis=0)
        RVstar_sorted = np.take_along_axis(self.RVstar, true_anomaly_sorted_idx, axis=0)
        vstar_sorted = np.take_along_axis(self.vstar, true_anomaly_sorted_idx, axis=0)  + RVstar_sorted
        vplanet_sorted = np.take_along_axis(self.RVplanet, true_anomaly_sorted_idx, axis=0)
        berv_sorted = np.take_along_axis(self.berv, true_anomaly_sorted_idx, axis=0)
        tel = - self.vsys + berv_sorted
        self.true_phases = self.true_anomaly_sorted / (2 * np.pi)
        self.transit_sorted =  np.take_along_axis(self.transit, true_anomaly_sorted_idx, axis=0)
        self.obstimes = np.take_along_axis(self.obstimes, true_anomaly_sorted_idx, axis=0)
    
        
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
            self.RVs_tel = tel

            print(f"[INFO] I don't know this restframe, I am using the system restframe")

