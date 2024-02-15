from RVTraceEstimator import RVTraceEstimator




W121 = RVTraceEstimator(dp='.', # datapath is just the same 
                        obs='paranal') 

W121.plot_doppler_trace(RF='star', save=True) # RF = rest frame


# access class elements

# e.g. RV_planet = W121.RVs_planet -- in the system above. So if you want it in a different system, change the restframe
# W121.transit gives you the tayph transit light curve and W121.true_phases gives you the orbital phase