from flask import Blueprint
from flask_restplus import Api

blueprint = Blueprint('main', __name__, url_prefix='/api')
api = Api(
    blueprint,
    title='Emission Simulation API',
    default='EM-VIZ BACKEND')

from app.main import routes