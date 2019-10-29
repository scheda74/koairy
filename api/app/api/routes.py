
from fastapi import APIRouter, Depends

from app.tools.simulation.parse_emission import Parser
from app.tools.simulation.simulator import Simulator
from app.tools.simulation.preprocessor import PreProcessor
from app.tools.regression.simple_lin_reg import LinReg
# from db.database import DB
from app.models.simulation_input import Inputs, example_body
# import db.query_database as query
from app.db.mongodb import AsyncIOMotorClient, get_database
from app.crud.emissions import get_caqi_emissions_for_sim

router = APIRouter()

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
    Starts a new simulation with given input parameters...
    """
    # body = await request.get_json()
    sim_id = generate_id(inputs)
    lr = LinReg(db, sim_id)
    await lr.predict_emission()
    # raw_emissions = await lr.fetch_simulated_emissions()
    # print(raw_emissions)
    # return raw_emissions["emissions"] if raw_emissions != None else {}
    return {}

def generate_id(inputs):
    src_weights = "".join([str(v).replace('.', '') for v in inputs.srcWeights.values()])
    dst_weights = "".join([str(v).replace('.', '') for v in inputs.dstWeights.values()])
    veh_dist = "".join([str(v).replace('.', '') for v in inputs.vehicleDistribution])
    return ("%s_%s_%s_%s_%s_%s" % (src_weights, dst_weights, veh_dist, inputs.vehicleNumber, inputs.timesteps, inputs.weatherScenario))
