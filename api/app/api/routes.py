
from fastapi import APIRouter

from tools.simulation.parse_emission import Parser
from tools.simulation.simulator import Simulator
from tools.simulation.preprocessor import PreProcessor
from db.database import DB
from models.simulation_input import Inputs, example_body
import db.query_database as query

router = APIRouter()

@router.post('/get/caqi')
async def get_caqi(inputs: Inputs = example_body):
    """
    Returns CAQI values. If not available new simulation will be started
    """
    
    # body = await request.get_json()
    simulation_id = generate_id(inputs)
    print(inputs)
    print(simulation_id)
    parser = Parser(simulation_id)
    # emission_cycle = await parser.parse_simulated_emissions()
    caqi = query.get_latest_emissions(simulation_id)
    if caqi != None:
        return caqi["emissions"] 
    else: 
        return parser.get_caqi_data()

# @router.get('/generate/weights')
# async def generate_weights(self):
#     """
#     Writes new weights with the given inputs from area distribution
#     """
#     body = await request.get_json()
#     simulation_id = generate_id(body)
#     processor = PreProcessor(
#         sim_id=simulation_id,
#         timesteps=body['timesteps'],
#         agents=body['vehicleNumber'], 
#         src_weights=body['srcWeights'], 
#         dst_weights=['dstWeights']
#     )
#     await processor.write_weight_file()
#     return "File written"

# @router.post('/start/simulation')
# async def start_simulation(self):
#     """
#     Starts a new simulation with given input parameters...
#     """
#     body = await request.get_json()
#     sim_id = generate_id(body)
#     print("Starting PreProcessor...")
#     processor = PreProcessor(
#         sim_id=sim_id,
#         timesteps=body['timesteps'],
#         agents=body['vehicleNumber'], 
#         src_weights=body['srcWeights'], 
#         dst_weights=body['dstWeights']
#     )
#     cfg_filepath = await processor.preprocess_simulation_input()
#     # print(cfg_filepath)
#     print("Starting SUMO...")
#     simulator = Simulator(cfg_filepath, sim_id)
#     await simulator.start()
    
#     print("Parsing results...")
#     parser = Parser(sim_id)
#     return parser.get_caqi_data()
    # return await parser.parse_emissions()
    # return "OK"

def generate_id(inputs):
    src_weights = "".join([str(v).replace('.', '') for v in inputs.srcWeights.values()])
    dst_weights = "".join([str(v).replace('.', '') for v in inputs.dstWeights.values()])
    veh_dist = "".join([str(v).replace('.', '') for v in inputs.vehicleDistribution])
    return ("%s_%s_%s_%s_%s_%s" % (src_weights, dst_weights, veh_dist, inputs.vehicleNumber, inputs.timesteps, inputs.weatherScenario))
