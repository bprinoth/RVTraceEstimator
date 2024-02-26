from flask import Flask, render_template, request, jsonify
from RVTraceEstimatorInteractive import RVTraceEstimator
import plotly.express as px
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
import numpy as np

app = Flask(__name__)

def retrieve_system_parameters(system_id, obs, RF):
    # Query the database and retrieve system parameters based on system ID
    parameters = NasaExoplanetArchive.query_criteria(table='pscomppars', where=f"pl_name='{str(system_id)}'")

    # Define a function to retrieve individual parameters or return default if not found
    def get_parameter(parameter_name, default_value):
        try:
            if len(parameters) == 1:
                try:
                    value = float(parameters[parameter_name][0].value)
                    if np.isnan(value):
                        return default_value
                    else: 
                        return value
                except:
                    value = float(parameters[parameter_name][0])
                    if np.isnan(value):
                        return default_value
                    else: 
                        return value
            else:
                raise ValueError("Multiple parameters found")
        except (KeyError, ValueError, TypeError):
            print(f"[INFO] Error occurred while retrieving parameter '{parameter_name}'. Defaulting to default value.")
            return default_value

    # Extract parameters or use defaults if not found
    transitC = get_parameter('pl_tranmid', default_value=0.0)
    period = get_parameter('pl_orbper', default_value=0.0)
    ecc = get_parameter('pl_orbeccen', default_value=0.0)
    omega = get_parameter('pl_orblper', default_value=0.0)
    a = get_parameter('pl_orbsmax', default_value=0.0)
    aRs = get_parameter('pl_ratdor', default_value=0.0)
    orbinc = get_parameter('pl_orbincl', default_value=0.0)
    vsini = get_parameter('st_vsin', default_value=0.0)
    pob = get_parameter('pl_projobliq', default_value=0.0)
    Ms = get_parameter('st_mass', default_value=0.0)
    Mp = get_parameter('pl_bmassj', default_value=0.0)
    RpRs = get_parameter('pl_ratror', default_value=0.0)
    K = get_parameter('pl_rvamp', default_value=0.0) / 1000
    vsys = get_parameter('st_radv', default_value=0.0)
    T14 = get_parameter('pl_trandur', default_value=0.0)
    RA = get_parameter('ra', default_value=0.0)
    DEC = get_parameter('dec', default_value=0.0)
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
    app.run(debug=True, port=5022)