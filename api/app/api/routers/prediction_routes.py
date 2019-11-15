
import os
import json
import datetime
import pandas as pd
import matplotlib.pyplot as plt
from fastapi import APIRouter, Depends
from starlette.responses import FileResponse, RedirectResponse

from app.tools.simulation.parse_emission import Parser
from app.tools.simulation.simulator import Simulator
from app.tools.simulation.preprocessor import SimulationPreProcessor
from app.tools.predictor.lin_reg import LinReg
from app.tools.predictor.neural_nets.py import NeuralNet
# from db.database import DB
from app.models.simulation_input import Inputs, example_body
from app.models.prediction_input import PlotInput, example_plot_input
# import db.query_database as query
from app.db.mongodb import AsyncIOMotorClient, get_database
from app.crud.emissions import get_caqi_emissions_for_sim
from app.crud.hawa_dawa import (
    get_hawa_dawa_by_time
    # get_all_hawa_dawa
)
from app.crud.bremicker import (
    get_bremicker
)
from app.core.config import PLOT_BASEDIR


router = APIRouter()

@router.post('/start/tbats')
async def start_tbats(inputs: Inputs = example_body, db: AsyncIOMotorClient=Depends(get_database)):
    sim_id = generate_id(inputs)
    lr = LinReg(db, sim_id)

    df_tbats = await lr.start_tbats(boxID=672, input_keys=['temp', 'hum', 'PMx', 'WIND_SPEED', 'WIND_DIR'], output_key='pm2.5')

# https://machinelearningmastery.com/time-series-prediction-lstm-recurrent-neural-networks-python-keras/
@router.post('/start/lstm')
async def start_lstm(inputs: Inputs = example_body, db: AsyncIOMotorClient=Depends(get_database)):
    sim_id = generate_id(inputs)
    nn = NeuralNet(db, sim_id)

    df_lstm = await nn.start_lstm(boxID=672, input_keys=['temp', 'hum', 'PMx', 'WIND_SPEED', 'WIND_DIR'], output_key='no2')


@router.post('/start/linreg')
async def start_linreg(inputs: Inputs = example_body, db: AsyncIOMotorClient=Depends(get_database)):
    """
    Starts training a simple Linear Regression Model with the specified independents (= inputs) and dependent (= output)
    Next, it'll predict the specified output with the data given to this request
    """
    sim_id = generate_id(inputs)
    lr = LinReg(db, sim_id)
    nn = NeuralNet(db, sim_id)
    df_lin = await lr.start_lin_reg(boxID=672, input_keys=['temp', 'hum', 'PMx', 'WIND_SPEED', 'WIND_DIR'], output_key='pm2.5')
    df_mlp = await nn.start_cnn(boxID=672, input_keys=['temp', 'hum', 'PMx', 'WIND_SPEED', 'WIND_DIR'], output_key='pm2.5')
    # df_no2_pred = await lr.predict_emission(boxID=672, input_keys=['temp', 'hum', 'NOx', 'WIND_SPEED', 'WIND_DIR'], output_key='no2')
    # print(df_pm10_pred)
    # print(df_no2_pred)
    # print(df_pm10_lin)
    df_combined = pd.concat([df_mlp, df_lin['pm2.5_lin_predicted']], axis=1)
    print(df_combined)
    # df_combined.plot(figsize=(18, 5))
    # plt.savefig(PLOT_BASEDIR + '/pm25_lin_mlp')
    # print(df_combined)
    # return df_combined.to_dict(orient='list')

    # return df_pm10_mlp

    # return await lr.get_sim_em_distribution()
    # return await lr.predict_emission()
    
    # result = await get_bremicker(db)
    # return result
    # result = await get_hawa_dawa_by_time(db)
    # print(result)
    # return result.to_dict(orient='list')
    # await lr.get_air_sensor_data()


@router.post('/start/cnn')
async def start_conv(inputs: Inputs = example_body, db: AsyncIOMotorClient=Depends(get_database)):
    sim_id = generate_id(inputs)
    nn = NeuralNet(db, sim_id)
    # res = await lr.aggregate_data(start_date='2019-08-01', end_date='2019-10-30', start_hour='0:00', end_hour='23:00')
    # return res.to_dict(orient='index')

    df_pm10_pred = await nn.start_cnn(boxID=672, input_keys=['temp', 'hum', 'PMx', 'WIND_SPEED', 'WIND_DIR'], output_key='pm10')
    print(df_pm10_pred)
    # df_pm25_pred = await lr.start_cnn(boxID=672, input_keys=['temp', 'hum', 'PMx'], output_key='pm2.5')
    # df_no2_pred = await lr.start_cnn(boxID=672, input_keys=['temp', 'hum', 'NOx'], output_key='no2')
    # df_combined = pd.concat([df_no2_pred, df_pm10_pred, df_pm25_pred], axis=1)
    # print(df_combined)
    # return df_combined.to_dict(orient='list')


def generate_id(inputs):
    src_weights = "".join([str(v).replace('.', '') for v in inputs.srcWeights.values()])
    dst_weights = "".join([str(v).replace('.', '') for v in inputs.dstWeights.values()])
    veh_dist = "".join([str(v).replace('.', '') for v in inputs.vehicleDistribution.values()])
    return ("%s_%s_%s_%s_%s_%s" % (src_weights, dst_weights, veh_dist, inputs.vehicleNumber, inputs.timesteps, inputs.weatherScenario))