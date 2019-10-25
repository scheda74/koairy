from quart_openapi import Pint, Resource, PintBlueprint, Swagger
from quart_cors import cors
from quart import request, jsonify
import json
# from simulation import parse_emission as parser
from simulation.parse_emission import Parser
from simulation.simulator import Simulator
from simulation.preprocessor import PreProcessor
from database.database import DB
import database.query_database as query

app = Pint(__name__, title="EM-ViZ API")

DB.init()
blueprint = PintBlueprint('main', __name__, url_prefix='/api')
blueprint = cors(blueprint, allow_origin="*")

expected = app.create_validator('simulation_request', {
    'weatherScenario': 'integer',
    'vehicleDistribution': [],
    'srcWeights': {'string': 'float'},
    'dstWeights': {'string': 'float'},
    'vehicleNumber': 'integer',
    'timesteps': 'integer'
    })

@blueprint.route('/')
class Root(Resource):
    async def get(self):
        '''
        Hello World Route

        This docstring will show up as the description and short-description
        for the openapi docs for this route.
        '''
        return "hello"

@blueprint.route('/get/emissions')
class Emission(Resource):
    async def get(self):
        """
        parses files and returns simulated emissions
        """
        # emission_cycle = await parser.parse_simulated_emissions()
        # emission_cycle = await parser.parse_emissions()
        # return emission_cycle
        return 'Emission'

@blueprint.route('/get/caqi', methods=['GET', 'POST'])
class CAQI(Resource):
    @blueprint.expect(expected)
    async def get(self):
        """
        Looks up database if same simulation has already been run
        If not: New simulation with input parameters will start
        """
        
        body = await request.get_json()
        simulation_id = generate_id(body)
        print(body)
        print(simulation_id)
        parser = Parser(simulation_id)
        # emission_cycle = await parser.parse_simulated_emissions()
        caqi = query.get_latest_emissions(simulation_id)
        if caqi != None:
            return caqi["emissions"] 
        else: 
            return parser.get_caqi_data()

@blueprint.route('/generate/weights')
class Weights(Resource):
    @blueprint.expect(expected)
    async def post(self):
        """
        Writes new weights with the given inputs from area distribution
        """
        body = await request.get_json()
        simulation_id = generate_id(body)
        processor = PreProcessor(
            sim_id=simulation_id,
            timesteps=body['timesteps'],
            agents=body['vehicleNumber'], 
            src_weights=body['srcWeights'], 
            dst_weights=['dstWeights']
        )
        await processor.write_weight_file()
        return "File written"

@blueprint.route('/start/simulation', methods=['GET', 'POST'])
class Simulation(Resource):
    @blueprint.expect(expected)
    async def get(self):
        """
        Starts a completely new simulation...
        """
        body = await request.get_json()
        sim_id = generate_id(body)
        print("Starting PreProcessor...")
        processor = PreProcessor(
            sim_id=sim_id,
            timesteps=body['timesteps'],
            agents=body['vehicleNumber'], 
            src_weights=body['srcWeights'], 
            dst_weights=body['dstWeights']
        )
        cfg_filepath = await processor.preprocess_simulation_input()
        # print(cfg_filepath)
        print("Starting SUMO...")
        simulator = Simulator(cfg_filepath, sim_id)
        await simulator.start()
        
        print("Parsing results...")
        parser = Parser(sim_id)
        return parser.get_caqi_data()
        # return await parser.parse_emissions()
        # return "OK"

def generate_id(data):
    src_weights = "".join([str(v).replace('.', '') for v in data['srcWeights'].values()])
    dst_weights = "".join([str(v).replace('.', '') for v in data['dstWeights'].values()])
    veh_dist = "".join([str(v).replace('.', '') for v in data['vehicleDistribution']])
    return ("%s_%s_%s_%s_%s_%s" % (src_weights, dst_weights, veh_dist, data['vehicleNumber'], data['timesteps'], data['weatherScenario']))

# Swagger(blueprint)
app.register_blueprint(blueprint)
app.run()