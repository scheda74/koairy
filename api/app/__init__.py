#!/usr/bin/env python
# coding=utf-8

# from flask import Flask
from quart import Quart
from quart_openapi import Pint

def create_app(config):
    # app = Flask(__name__)
    # app = Pint(__name__, title="EM-VIZ BACKEND")
    app = Quart(__name__)
    # register_blueprints(app)
    return app


# def register_blueprints(app):
#     from app.main import blueprint as main_api
#     app.register_blueprint(main_api)