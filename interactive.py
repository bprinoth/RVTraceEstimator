from flask import Flask, render_template, request, jsonify, send_file
from RVTraceEstimatorInteractive import RVTraceEstimator
import plotly.express as px
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
import numpy as np
from astropy.time import Time
from datetime import datetime, timezone
import pandas as pd


app = Flask(__name__)

def retrieve_system_parameters(system_id, obs, RF, obsdate, T14):
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
    observation_date = Time(obsdate, scale='utc').isot
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
    #T14 = get_parameter('pl_trandur', default_value=0.0)
    RA = get_parameter('ra', default_value=0.0)
    DEC = get_parameter('dec', default_value=0.0)
    n_exp = 100

    # Return the extracted values as a dictionary
    return {
        'observation_date': observation_date,
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
    obsdate = request.form['observation_date']
    T14 = request.form['T14']
    system_parameters = retrieve_system_parameters(system_id, obs, RF, obsdate, T14)
    return jsonify(system_parameters)

@app.route('/plot', methods=['POST'])
def plot():
    input_values = request.json['input_values']
    obs = input_values['obs']
    RF = input_values['RF']
    obsdate = input_values['observation_date']
    T14 = input_values['T14']

    input_values = {key: value for key, value in input_values.items() if key not in ['obs', 'RF', 'observation_date', 'T14']}

    input_params = []
    for key, value in input_values.items():
        if isinstance(value, str):
            input_params.append(float(value))

    rv = RVTraceEstimator(input_params=input_params, obs=obs, obsdate=obsdate, Tdur=T14)
    rv.calculate_RV_traces(RF=RF)

    y_axis = rv.true_phases

    fig = px.scatter()
    fig.add_scatter(x=list(rv.RVs_RM[rv.transit_sorted != 1]), y=list(y_axis[rv.transit_sorted != 1]), mode='lines', name='RM')
    fig.add_scatter(x=list(rv.RVs_planet.value), y=list(y_axis), mode='lines', name='Planet')
    fig.add_scatter(x=list(rv.RVs_star), y=list(y_axis), mode='lines', name='Star')
    fig.add_scatter(x=list(rv.RVs_tel), y=list(y_axis), mode='lines', name='Tellurics')

    fig.update_layout(
        xaxis_title="RVs [km/s]",
        yaxis_title=fr"Orbital phase"
    )

    phases = rv.true_phases.tolist()
    RVs_RM = rv.RVs_RM.tolist()
    transit_sorted = rv.transit_sorted.tolist()
    RVs_planet = rv.RVs_planet.value.tolist()
    RVs_star = rv.RVs_star.tolist()
    RVs_tel = rv.RVs_tel.tolist()

    return jsonify({
        'plot': fig.to_dict(),
        'phases': phases,
        'RVs_RM': RVs_RM,
        'transit_sorted': transit_sorted,
        'RVs_planet': RVs_planet,
        'RVs_star': RVs_star,
        'RVs_tel': RVs_tel
    })


@app.route('/download', methods=['POST'])

def download():
    data = request.json
    phases = data['phases']
    RVs_RM = data['RVs_RM']
    transit_sorted = data['transit_sorted']
    RVs_planet = data['RVs_planet']
    RVs_star = data['RVs_star']
    RVs_tel = data['RVs_tel']

    # Create a DataFrame
    df = pd.DataFrame({
        'phases': phases,
        'RVs_RM': RVs_RM,
        'transit_sorted': transit_sorted,
        'RVs_planet': RVs_planet,
        'RVs_star': RVs_star,
        'RVs_tel': RVs_tel
    })

    # Save the DataFrame to a CSV file
    df.to_csv('output.csv', index=False)

    # Return the file as an attachment
    return send_file('output.csv', as_attachment=True)


if __name__ == '__main__':
    app.run(debug=True, port=5021)