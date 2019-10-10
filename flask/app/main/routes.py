from flask import jsonify
from flask_restplus import Resource
from flask import request
from flask import make_response

from app.main import api
import docker
import datetime
import io
import requests
import logging

logger = logging.getLogger(__name__)

@api.route('/client/traffic/historic')
class HistoricTraffic(Resource):
    @staticmethod
    def get():
        # TODO: implement method; query database for historic data
        """
        returns historic traffic data
        """
        return 'Historic traffic'


# This route gets real-time traffic data from MongoDB and passes it on
@api.route('/client/traffic/realtime')
class RealtimeTraffic(Resource):
    @staticmethod
    def get():
        """
        returns realtime traffic data
        """
        return 'real time traffic'


@api.blueprint.after_request
def after_request(response):
    header = response.headers
    header['Access-Control-Allow-Origin'] = '*'
    return response