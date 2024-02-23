import dash 
from dash import dcc
from dash import html
from dash import Dash
from dash import State
from dash.dependencies import Input, Output, State
import plotly.express as px
import pandas as pd
import astropy.units as u
import numpy as np
from RVTraceEstimatorInteractive import RVTraceEstimator
import plotly.graph_objects as go
import dash_bootstrap_components as dbc
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive


def GENERATE_TRACES():
    app = dash.Dash(__name__, external_stylesheets=[dbc.themes.FLATLY])
    
    # Function to retrieve system parameters based on system ID
    def retrieve_system_parameters(system_id):
        #print(system_id)
        # Your code to query the database and retrieve system parameters based on system ID
        # For demonstration purposes, I'll just return some default values
        parameters = NasaExoplanetArchive.query_criteria(table='pscomppars', where=f"pl_name='{str(system_id)}'")
    
        #print(len(parameters))
    
        if len(parameters) == 1:
            
            #parameters = planet_data[0]
            transitC = parameters['pl_tranmid']
            
            period = parameters['pl_orbper']
            ecc = parameters['pl_orbeccen']
            omega = parameters['pl_orblper']
            a = parameters['pl_orbsmax']
            aRs = parameters['pl_ratdor']
            orbinc = parameters['pl_orbincl']
            vsini = parameters['st_vsin']
            pob = parameters['pl_projobliq']
            Ms = parameters['st_mass']
            Mp = parameters['pl_bmassj']
            RpRs = parameters['pl_ratror']
            K = parameters['pl_rvamp'] / 1000
            vsys = parameters['st_radv']
            T14 = parameters['pl_trandur']
            RA = parameters['ra']
            DEC = parameters['dec']
            n_exp = 100  # Assuming a default value
            #return (transitC[0].value, period, ecc, omega, a, aRs, orbinc, vsini, pob, Ms, Mp, RpRs, K, vsys, T14, RA, 
    
            return {'transitC': transitC[0].value,
                'period': period[0].value,
                'ecc': ecc[0],
                'omega': omega[0].value,
                'a': a[0].value,
                'aRs': aRs[0],
                'orbinc': orbinc[0].value,
                'vsini': vsini[0].value,
                'pob': pob[0].value,
                'Ms': Ms[0].value,
                'Mp': Mp[0].value,
                'RpRs': RpRs[0],
                'K': K[0].value,
                'vsys': vsys[0].value,
                'T14': T14[0].value,
                'RA': RA[0].value,
                'DEC': DEC[0].value,
                'n_exp': 100
                
            }
    
        else:
            return {
                'transitC': 0,
                'period': 0,
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
                'n_exp': 100
            }
    
    app.layout = html.Div([
        dbc.Row([
            dbc.Col([
                dcc.Input(placeholder='System ID', type='text', id='system-id-input'),
                html.Button('Retrieve Parameters', id='retrieve-parameters-button', className='btn btn-primary'),
                html.Button('Plot', id='plot-button', className='btn btn-success')
            ], width=6)
        ]),
        dbc.Row([
            dbc.Col([
                html.Div(id='parameter-inputs', children=[
                    html.Div([
                        html.Label('Transit centre time in BJD of your observation night'),
                        dcc.Input(placeholder='', type='number', id='transitC-value')
                    ]),
                    html.Div([
                        html.Label('Period in days'),
                        dcc.Input(placeholder='', type='number', id='period-value')
                    ]),
                    html.Div([
                        html.Label('Eccentricity'),
                        dcc.Input(placeholder='', type='number', id='ecc-value')
                    ]),
                    html.Div([
                        html.Label(r'Argument of periastron $\omega$ in degrees'),
                        dcc.Input(placeholder='', type='number', id='omega-value')
                    ]),
                    html.Div([
                        html.Label('Semi-major axis in AU'),
                        dcc.Input(placeholder='', type='number', id='a-value')
                    ]),
                    html.Div([
                        html.Label('Scaled semi-major axis'),
                        dcc.Input(placeholder='', type='number', id='aRs-value')
                    ]),
                    html.Div([
                        html.Label('Orbital inclination in degrees'),
                        dcc.Input(placeholder='', type='number', id='orbinc-value')
                    ]),
                    html.Div([
                        html.Label('vsini (km/s)'),
                        dcc.Input(placeholder='', type='number', id='vsini-value')
                    ]),
                    html.Div([
                        html.Label(r'Projected orbital obliquity $\lambda$'),
                        dcc.Input(placeholder='', type='number', id='pob-value')
                    ]),
                    html.Div([
                        html.Label('Stellar mass [Msun]'),
                        dcc.Input(placeholder='', type='number', id='Ms-value')
                    ]),
                    html.Div([
                        html.Label('Planetary mass [Mjup]'),
                        dcc.Input(placeholder='', type='number', id='Mp-value')
                    ]),
                ]),
            ]),
            dbc.Col([
                html.Div(id='parameter-inputs2', children=[
                    html.Div([
                        html.Label('Scaled planetary radius'),
                        dcc.Input(placeholder='', type='number', id='RpRs-value')
                    ]),
                    html.Div([
                        html.Label('Radial velocity amplitude K (m/s)'),
                        dcc.Input(placeholder='', type='number', id='K-value')
                    ]),
                    html.Div([
                        html.Label('Systemic velocity vsys (m/s)'),
                        dcc.Input(placeholder='', type='number', id='vsys-value')
                    ]),
                    html.Div([
                        html.Label('Transit duration T14 (hours)'),
                        dcc.Input(placeholder='', type='number', id='T14-value')
                    ]),
                    html.Div([
                        html.Label('Right Ascension RA (degrees)'),
                        dcc.Input(placeholder='', type='text', id='RA-value')
                    ]),
                    html.Div([
                        html.Label('Declination DEC (degrees)'),
                        dcc.Input(placeholder='', type='text', id='DEC-value')
                    ]),
                    html.Div([
                        html.Label('Number of exposures per orbit'),
                        dcc.Input(placeholder='', type='number', id='n_exp-value')
                    ]),
                    dcc.Dropdown(
                        options=[{'label': 'paranal', 'value': 'paranal'}], 
                        value='paranal',
                        multi=False,
                        id='obs-drop'
                    ),
                    dcc.Dropdown(
                        options=[
                            {'label': 'system', 'value': 'system'},
                            {'label': 'star', 'value': 'star'},
                            {'label': 'planet', 'value': 'planet'},
                            {'label': 'obs', 'value': 'obs'},
                            {'label': 'berv', 'value': 'berv'}
                        ],
                        value='obs',
                        multi=False,
                        id='RF-drop',
                    )
                ])
            ])
        ]),
        dcc.Graph(id='indicator-graphic')
    ])
    
    
    @app.callback(
        Output('parameter-inputs', 'children'),
        [Input('retrieve-parameters-button', 'n_clicks')],
        [State('system-id-input', 'value')]
    )
    def retrieve_parameters(n_clicks, system_id):
        if n_clicks:
            system_parameters = retrieve_system_parameters(system_id)
        else:
            system_parameters = {}
    
        parameter_inputs = [
            html.Div([html.Label('Transit centre time in BJD of your observation night')]),
            dcc.Input(placeholder='', type='number', id='transitC-value', value=system_parameters.get('transitC', '')),
            html.Div([html.Label('Period in days')]),
            dcc.Input(placeholder='', type='number', id='period-value', value=system_parameters.get('period', '')),
            html.Div([html.Label('Eccentricity')]),
            dcc.Input(placeholder='', type='number', id='ecc-value', value=system_parameters.get('ecc', '')),
            html.Div([html.Label(r'Argument of periastron $\omega$ in degrees')]),
            dcc.Input(placeholder='', type='number', id='omega-value', value=system_parameters.get('omega', '')),
            html.Div([html.Label('Semi-major axis in AU')]),
            dcc.Input(placeholder='', type='number', id='a-value', value=system_parameters.get('a', '')),
            html.Div([html.Label('Scaled semi-major axis')]),
            dcc.Input(placeholder='', type='number', id='aRs-value', value=system_parameters.get('aRs', '')),
            html.Div([html.Label('Orbital inclination in degrees')]),
            dcc.Input(placeholder='', type='number', id='orbinc-value', value=system_parameters.get('orbinc', '')),
            html.Div([html.Label('vsini (km/s)')]),
            dcc.Input(placeholder='', type='number', id='vsini-value', value=system_parameters.get('vsini', '')),
            html.Div([html.Label(r'Projected orbital obliquity $\lambda$')]),
            dcc.Input(placeholder='', type='number', id='pob-value', value=system_parameters.get('pob', '')),
            html.Div([html.Label('Stellar mass [Msun]')]),
            dcc.Input(placeholder='', type='number', id='Ms-value', value=system_parameters.get('Ms', '')),
            html.Div([html.Label('Planetary mass [Mjup]')]),
            dcc.Input(placeholder='', type='number', id='Mp-value', value=system_parameters.get('Mp', ''))
        ]
        return parameter_inputs
    
    @app.callback(
        Output('parameter-inputs2', 'children'),
        [Input('retrieve-parameters-button', 'n_clicks')],
        [State('system-id-input', 'value')]
    )
    def retrieve_parameters2(n_clicks, system_id):
        if n_clicks:
            system_parameters = retrieve_system_parameters(system_id)
        else:
            system_parameters = {}
    
        parameter_inputs = [
            html.Div([html.Label('Scaled planetary radius')]),
            dcc.Input(placeholder='', type='number', id='RpRs-value', value=system_parameters.get('RpRs', '')),
            html.Div([html.Label('Radial velocity amplitude K (m/s)')]),
            dcc.Input(placeholder='', type='number', id='K-value', value=system_parameters.get('K', '')),
            html.Div([html.Label('Systemic velocity vsys (m/s)')]),
            dcc.Input(placeholder='', type='number', id='vsys-value', value=system_parameters.get('vsys', '')),
            html.Div([html.Label('Transit duration T14 (hours)')]),
            dcc.Input(placeholder='', type='number', id='T14-value', value=system_parameters.get('T14', '')),
            html.Div([html.Label('Right Ascension RA (degrees)')]),
            dcc.Input(placeholder='', type='text', id='RA-value', value=system_parameters.get('RA', '')),
            html.Div([html.Label('Declination DEC (degrees)')]),
            dcc.Input(placeholder='', type='text', id='DEC-value', value=system_parameters.get('DEC', '')),
            html.Div([html.Label('Number of exposures per orbit')]),
            dcc.Input(placeholder='', type='number', id='n_exp-value', value=system_parameters.get('n_exp', '')),
            html.Div([html.Label('Observatory')]),
                dcc.Dropdown(
                    options=[{'label': 'paranal', 'value': 'paranal'}], 
                    value='paranal',
                    multi=False,
                    id='obs-drop'
                ),
                html.Div([html.Label('Reference Frame')]),
                dcc.Dropdown(
                    options=[
                        {'label': 'system', 'value': 'system'},
                        {'label': 'star', 'value': 'star'},
                        {'label': 'planet', 'value': 'planet'},
                        {'label': 'obs', 'value': 'obs'},
                        {'label': 'berv', 'value': 'berv'}
                    ],
                    value='obs',
                    multi=False,
                    id='RF-drop',
                )
            ]
        return parameter_inputs
    
    
    # Callback to update the graph based on input parameters
    @app.callback(
        Output('indicator-graphic', 'figure'),
        [Input('plot-button', 'n_clicks')],
        [State('parameter-inputs', 'children'),
         State('parameter-inputs2', 'children')]
    )
    def update_graph(n_clicks, parameter_inputs1, parameter_inputs2):
        if n_clicks:
            input_values = []
            all_parameter_inputs = parameter_inputs1 + parameter_inputs2
            for item in all_parameter_inputs:
                if item['type'] == 'Input':
                    input_values.append(item['props']['value'])
                elif item['type'] == 'Dropdown':
                    input_values.append(item['props']['value'])
    
            #print(input_values)
            
            rv = RVTraceEstimator(input_params=input_values[:-2:], obs=input_values[-2])  # Assuming 'obs' is default observation
            rv.calculate_RV_traces(RF=input_values[-1])  # Assuming RF is already defined somewhere
            
            #colors = sns.color_palette('colorblind', 4)
    
            fig = fig = px.scatter()
            fig.add_scatter(x=rv.RVs_RM[rv.transit_sorted != 1], y=rv.true_phases[rv.transit_sorted != 1], mode='lines', name='RM')
            fig.add_scatter(x=rv.RVs_planet, y=rv.true_phases, mode='lines', name='Planet')
            fig.add_scatter(x=rv.RVs_star, y=rv.true_phases, mode='lines', name='Star')
            fig.add_scatter(x=rv.RVs_tel, y=rv.true_phases, mode='lines', name='Tellurics')
            
            return fig
        else:
            return dash.no_update

    return app