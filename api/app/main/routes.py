import asyncio
import json
from quart import Quart, request, url_for, jsonify, Blueprint
from quart_cors import cors
from quart_openapi import Resource
# from app.simulation import parse_emission as parser
from app.simulation import parse_emission as parser
from app.simulation.simulator import Simulator
from app.simulation import preprocessor as ip
from app.database.database import DB
import app.database.query_database as query
from app.main import blueprint

@blueprint.route('/')
class Root(Resource):
    async def get(self):
        '''Hello World Route

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

@blueprint.route('/get/caqi')
class CAQI(Resource):
    async def get(self):
        """
        Looks up database if same simulation has already been run
        If not: New simulation with input parameters will start
        """
        # emission_cycle = await parser.parse_simulated_emissions()
        caqi = query.get_latest_emissions()
        if caqi != None:
            return caqi["emissions"] 
        else: 
            return parser.get_caqi_data()

@blueprint.route('/generate/weights')
class Weights(Resource):
    async def post(self):
        """
        Writes new weights with the given inputs from area distribution
        """
        processor = ip.PreProcessor()
        await processor.write_weight_file()
        return "File written"

@blueprint.route('/start/simulation')
class Simulation(Resource):
    async def get(self):
        """
        Starts a completely new simulation...
        """
        print("Starting PreProcessor...")
        processor = ip.PreProcessor()
        cfg_filepath = await processor.preprocess_simulation_input()

        print("Starting SUMO...")
        simulator = Simulator(cfg_filepath)
        await simulator.start()
        
        print("Parsing results...")
        return parser.get_caqi_data()
        # return await parser.parse_emissions()
        # return "OK"

# if __name__ == "__main__":
#     app.run('localhost', port=5000, debug=True)