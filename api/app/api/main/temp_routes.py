# from flask import jsonify
from flask_restplus import Resource
from flask import request
from flask import make_response

# from quart import request
# from quart import make_response
# from quart_openapi import Pint, Resource

from app.main import api
from app.simulation import parse_emission as parser
from app.simulation import simulator
from app.main import blueprint

# import docker
# import datetime
# import io
# import requests
# import logging
import asyncio

# logger = logging.getLogger(__name__)

@api.route('/get/emissions')
# @app.route('/get/emissions')
# @blueprint.route('/get/emissions')
class SimulatedEmissions(Resource):
    @staticmethod
    def get():
        # TODO: implement method; query database for historic data
        """
        parses files and returns simulated emissions
        """
        # emission_cycle = await parser.parse_emissions()
        # return emission_cycle
        return 'Emissions'
    # blueprint = Blueprint(__name__)
    # blueprint.add_url_rule("/get/emissions", get)

# This route gets real-time traffic data from MongoDB and passes it on
# @blueprint.route('/client/traffic/realtime')
@api.route('/client/traffic/realtime')
class RealtimeTraffic(Resource):
    @staticmethod
    def get():
        """
        returns realtime traffic data
        """
        return 'real time traffic'

@api.route('/get/caqi')
class AirQuality(Resource):
    @staticmethod
    def get():
        """
        parses files and returns simulated emissions as CAQI
        """
        # emission_cycle = await parser.parse_simulated_emissions()
        return parser.get_caqi_data()

@api.route('/start/simulation')
class Simulation(Resource):
    @staticmethod
    async def get():
        """
        parses files and returns simulated emissions
        """
        # simulator.start()
        print("now I'm parsing results")
        # await simulator.start()
        asyncio.gather(*simulator.start())
        print("Now I'm using the emissions from simulation")
        result = asyncio.gather(*parser.get_caqi_data())
        # return await parser.get_caqi_data()
        return result
        # return await parser.parse_simulated_emissions()
        # return parser.parse_emissions()
        # return "OK"

@api.blueprint.after_request
def after_request(response):
    header = response.headers
    header['Access-Control-Allow-Origin'] = '*'
    return response