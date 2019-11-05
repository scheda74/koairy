
import os
import json
import datetime
import pandas as pd
import matplotlib.pyplot as plt
from fastapi import APIRouter, Depends
from starlette.responses import FileResponse, RedirectResponse

from app.tools.simulation.parse_emission import Parser
from app.tools.simulation.simulator import Simulator
from app.tools.simulation.preprocessor import PreProcessor
from app.tools.regression.simple_lin_reg import LinReg
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

@router.get("/")
async def redirect():
    response = RedirectResponse(url='/docs')
    return response

@router.post('/get/caqi')
async def get_caqi(inputs: Inputs = example_body, db: AsyncIOMotorClient=Depends(get_database)):
    """
    Returns CAQI values. If not available new simulation will be started
    """
    print(inputs)
    simulation_id = generate_id(inputs)
    print("[PARSER] Get CAQI data from simulation with id {simulation_id}")
    parser = Parser(db, simulation_id)
    return await parser.get_caqi_data()

@router.get('/generate/weights')
async def generate_weights(inputs: Inputs = example_body, db: AsyncIOMotorClient=Depends(get_database)):
    """
    Writes new weights with the given inputs from area distribution
    """
    sim_id = generate_id(inputs)
    processor = PreProcessor(
        sim_id=sim_id,
        timesteps=inputs.timesteps,
        agents=inputs.vehicleNumber, 
        src_weights=inputs.srcWeights, 
        dst_weights=inputs.dstWeights
    )
    await processor.write_weight_file()
    return "File written"

@router.post('/start/simulation')
async def start_simulation(inputs: Inputs = example_body, db: AsyncIOMotorClient=Depends(get_database)):
    """
    Starts a new simulation with given input parameters...
    """
    # body = await request.get_json()
    sim_id = generate_id(inputs)
    print("Starting PreProcessor...")
    processor = PreProcessor(
        sim_id=sim_id,
        timesteps=inputs.timesteps,
        agents=inputs.vehicleNumber, 
        src_weights=inputs.srcWeights, 
        dst_weights=inputs.dstWeights
    )
    cfg_filepath = await processor.preprocess_simulation_input()
    # print(cfg_filepath)
    print("Starting SUMO...")
    simulator = Simulator(db, cfg_filepath, sim_id)
    await simulator.start()
    
    print("Parsing results...")
    parser = Parser(db, sim_id)
    return await parser.get_caqi_data()
    # return await parser.parse_emissions()
    # return "OK"

@router.post('/start/linreg')
async def start_linreg(inputs: Inputs = example_body, db: AsyncIOMotorClient=Depends(get_database)):
    """
    Starts training a simple Linear Regression Model with the specified independents (= inputs) and dependent (= output)
    Next, it'll predict the specified output with the data given to this request
    """
    sim_id = generate_id(inputs)
    lr = LinReg(db, sim_id, boxID=672)
    df_pm10_pred = await lr.predict_emission(boxID=672, input_keys=['temp', 'hum', 'PMx'], output_key='pm10')
    df_no2_pred = await lr.predict_emission(boxID=672, input_keys=['temp', 'hum', 'NOx'], output_key='no2')
    print(df_pm10_pred)
    print(df_no2_pred)
    df_combined = pd.concat([df_no2_pred, df_pm10_pred], axis=1)
    print(df_combined)
    return df_combined.to_dict(orient='list')

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
    lr = LinReg(db, sim_id, boxID=672)
    await lr.start_cnn()

@router.post('/get/mean/vehicle')
async def get_mean_vehicles(inputs: Inputs = example_body, db: AsyncIOMotorClient=Depends(get_database)):
    """
    Starts training a simple Linear Regression Model with the specified independents (= inputs) and dependent (= output)
    Next, it'll predict the specified output with the data given to this request
    """
    sim_id = generate_id(inputs)
    lr = LinReg(db, sim_id)
    df = await lr.get_mean_vehicle_by_hour('2019-08-22', '2019-10-24', '07:00', '11:00')
    return df.to_dict(orient='list')

@router.post('/get/sensors')
async def get_sensors(db: AsyncIOMotorClient=Depends(get_database)):
    lr = LinReg(db)
    # await lr.get_hw_data()

@router.post('/get/plot')
async def get_plot(inputs: PlotInput = example_plot_input, db: AsyncIOMotorClient=Depends(get_database)):
    """
    Plot a line chart of the given attributes in a specific timeframe
    """
    lr = LinReg(db)
    df = await lr.aggregate_data(inputs.start_date, inputs.end_date, inputs.start_hour, inputs.end_hour)
    df = df[inputs.keys_to_compare]

    filename = PLOT_BASEDIR + '/' + 'fromrequest'
    df.plot(figsize=(18, 5))
    plt.ioff()
    plt.savefig(filename)
    # return FileResponse(plt.savefig(), media_type='image/png')
    return 'Done'


def generate_id(inputs):
    src_weights = "".join([str(v).replace('.', '') for v in inputs.srcWeights.values()])
    dst_weights = "".join([str(v).replace('.', '') for v in inputs.dstWeights.values()])
    veh_dist = "".join([str(v).replace('.', '') for v in inputs.vehicleDistribution])
    return ("%s_%s_%s_%s_%s_%s" % (src_weights, dst_weights, veh_dist, inputs.vehicleNumber, inputs.timesteps, inputs.weatherScenario))
