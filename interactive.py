from flask import Flask, render_template, request, jsonify
from RVTraceEstimatorInteractive import RVTraceEstimator
import plotly.express as px
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
import numpy as np

app = Flask(__name__)


def retrieve_system_parameters(system_id, obs, RF):
    # Query the database and retrieve system parameters based on system ID
    parameters = NasaExoplanetArchive.query_criteria(table='pscomppars', where=f"pl_name='{str(system_id)}'")
    
    if len(parameters) == 1:
        # Extract parameter values
        transitC = float(parameters['pl_tranmid'][0].value)
        period = float(parameters['pl_orbper'][0].value) # if 'pl_orbper' in parameters else np.nan
        ecc = float(parameters['pl_orbeccen'][0]) # if 'pl_orbeccen' in parameters else np.nan
        omega = float(parameters['pl_orblper'][0].value)# if 'pl_orblper' in parameters else np.nan
        a = float(parameters['pl_orbsmax'][0].value)# if 'pl_orbsmax' in parameters else np.nan
        aRs = float(parameters['pl_ratdor'][0]) # if 'pl_ratdor' in parameters else np.nan
        orbinc = float(parameters['pl_orbincl'][0].value) # if 'pl_orbincl' in parameters else np.nan
        vsini = float(parameters['st_vsin'][0].value) #if 'st_vsin' in parameters else np.nan
        pob = float(parameters['pl_projobliq'][0].value) #if 'pl_projobliq' in parameters else np.nan
        Ms = float(parameters['st_mass'][0].value) #if 'st_mass' in parameters else np.nan
        Mp = float(parameters['pl_bmassj'][0].value) #if 'pl_bmassj' in parameters else np.nan
        RpRs = float(parameters['pl_ratror'][0]) #if 'pl_ratror' in parameters else np.nan
        K = float(parameters['pl_rvamp'][0].value / 1000) #if 'pl_rvamp' in parameters else np.nan
        vsys = float(parameters['st_radv'][0].value) #if 'st_radv' in parameters else np.nan
        T14 = float(parameters['pl_trandur'][0].value) #if 'pl_trandur' in parameters else np.nan
        RA = float(parameters['ra'][0].value) # if 'ra' in parameters else np.nan
        DEC = float(parameters['dec'][0].value) #if 'dec' in parameters else np.nan
        n_exp = 100
        
        # Return the extracted values as a dictionary
        return {
            'transitC': transitC,
            'period': period,
            'ecc': ecc,
            'omega': omega,
            'a': a,
            'aRs': aRs,
            'orbinc': orbinc,
            'vsini': vsini,
            'pob': pob,
            'Ms': Ms,
            'Mp': Mp,
            'RpRs': RpRs,
            'K': K,
            'vsys': vsys,
            'T14': T14,
            'RA': RA,
            'DEC': DEC,
            'n_exp': n_exp,
            'obs': obs,
            'RF': RF
        }
    else:
        # Return default values if no parameters are found
        return {
            'transitC': np.nan,
            'period': np.nan,
            'ecc': 0.7,
            'omega': 60,
            'a': 0.155,
            'aRs': 18,
            'orbinc': 90,
            'vsini': 20,
            'pob': 1.2,
            'Ms': 1.53,
            'Mp': 4.0,
            'RpRs': 0.07,
            'K': 0.338,
            'vsys': 20,
            'T14': 18,
            'RA': '10:23:56.2397428368',
            'DEC': '-56:50:35.342013012',
            'n_exp': 100,
            'obs': 'paranal',
            'RF': 'star'
        }

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/get_parameters', methods=['POST'])
def get_parameters():
    system_id = request.form['system_id']
    obs = request.form['obs']
    RF = request.form['RF']
    system_parameters = retrieve_system_parameters(system_id, obs, RF)
    return jsonify(system_parameters)

@app.route('/plot', methods=['POST'])
def plot():
    input_values = request.json['input_values']
    #print(input_values)
    # Extract observatory and rest frame without converting to strings
    obs = input_values['obs'] 
    RF = input_values['RF']

    #print(obs, RF)
    # Remove observatory and rest frame from input values
    input_params = [value for key, value in input_values.items() if key not in ['obs', 'RF']]

    # Convert numeric values to float, leave non-numeric values unchanged
    input_params = [float(value.replace(',', '')) if isinstance(value, str) and ',' in value else float(value) for value in input_params]

    rv = RVTraceEstimator(input_params=input_params, obs=obs)
    rv.calculate_RV_traces(RF=RF)
    
    y_axis = (rv.obstimes - rv.Tc_observation) * 24 # to make it into hours
    
    centre_phase = rv.true_phases[np.argmin(np.abs(y_axis))]
    
    #print(y_axis)

    fig = px.scatter()
    fig.add_scatter(x=rv.RVs_RM[rv.transit_sorted != 1], y=y_axis[rv.transit_sorted != 1], mode='lines', name='RM')
    fig.add_scatter(x=rv.RVs_planet, y=y_axis, mode='lines', name='Planet')
    fig.add_scatter(x=rv.RVs_star, y=y_axis, mode='lines', name='Star')
    fig.add_scatter(x=rv.RVs_tel, y=y_axis, mode='lines', name='Tellurics')
    
    fig.update_layout(
        xaxis_title="RVs [km/s]",
        yaxis_title=fr"Time from centre phase $\phi$ = {np.round(centre_phase,2)} [h]"
    )

    return fig.to_json()



if __name__ == '__main__':
    app.run(debug=True, port=5021)