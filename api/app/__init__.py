# # # coding=utf-8

# # from flask import Flask
# from quart import Quart
# from quart_cors import cors
# from quart_openapi import Pint
# from app.database.database import DB
# # # from quart_openapi import Pint

# # def create_app(config):
# #     # app = Flask(__name__)
# #     # app = Pint(__name__, title="EM-VIZ BACKEND")
# #     app = Flask(__name__)
# #     register_blueprints(app)
# #     return app


# # def register_blueprints(app):
# #     from app.main import blueprint as main_api
# #     app.register_blueprint(main_api)

# # app = Quart(__name__)
# # app = cors(app, allow_origin="*")
# # DB.init()

# def create_app(config):
#     app = Pint(__name__, title="EM-ViZ API")
#     DB.init()
#     register_blueprints(app)
#     app.run()
#     return app

# # from app.main import routes
# # from app.main import blueprint as main_api
# # app.register_blueprint(main_api)

# def register_blueprints(app):
#     from app.main import blueprint as main_api
#     app.register_blueprint(main_api)

