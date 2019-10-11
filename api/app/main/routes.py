# from flask import jsonify
# from flask_restplus import Resource
# from flask import request
# from flask import make_response

# from quart import request
# from quart import make_response
from quart_openapi import Pint, Resource

# from app.main import api
from app.simulation import parse_emission as parser
from app.main import blueprint

# import docker
# import datetime
# import io
# import requests
# import logging
# import asyncio

# logger = logging.getLogger(__name__)

# @api.route('/get/emissions')
# @app.route('/get/emissions')
@blueprint.route('/get/emissions')
class SimulatedEmissions(Resource):
    @staticmethod
    async def get():
        # TODO: implement method; query database for historic data
        """
        parses files and returns simulated emissions
        """
        emission_cycle = await parser.parse_simulated_emissions()
        return emission_cycle
        # return 'Historic traffic'
    # blueprint = Blueprint(__name__)
    # blueprint.add_url_rule("/get/emissions", get)

# This route gets real-time traffic data from MongoDB and passes it on
# @app.route('/client/traffic/realtime')
@blueprint.route('/client/traffic/realtime')
class RealtimeTraffic(Resource):
    @staticmethod
    def get():
        """
        returns realtime traffic data
        """
        return 'real time traffic'


@blueprint.after_request
def after_request(response):
    header = response.headers
    header['Access-Control-Allow-Origin'] = '*'
    return response