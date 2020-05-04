"""This is the main file that runs the app. FLASK_APP env variable should be
set to a path to this file (e.g. on Windows "set FLASK_APP=api/run.py")."""

import logging
import os
from logging.handlers import RotatingFileHandler

from flask import Flask
from flask.logging import default_handler
from flask_cors import CORS

import sys


sys.path.insert(0, 'api')  # required in deployment to import api modules
sys.path.insert(0, '../../MINE-Database')
sys.path.insert(0, '..')  # required in deployment to import api modules


from api.config import Config
from api.database import mongo
from api.routes import mineserver_api


def create_app(instance_config=Config):
    """Create a Flask instance of MINE-Server with a specified configuration.

    Parameters
    ----------
    instance_config : Config (default: app.config.Config)
        Specifies configuration for this instance. Useful to change when
        running unit tests (see app.config_unittest.ConfigUnitTest).

    Notes
    -----
    This method is automatically called by Flask upon startup.
    """

    # Initialize app
    app = Flask(__name__)
    app.config.from_object(instance_config)

    # Register routes
    app.register_blueprint(mineserver_api, url_prefix='/mineserver')

    # Connect to Mongo Database
    mongo.init_app(app)

    # Allow CORS so we can have front end and back end on same server
    CORS(app)

    # Initialize logger
    if __name__ != '__main__':
        gunicorn_logger = logging.getLogger('gunicorn.error')
        app.logger.handlers = gunicorn_logger.handlers
        app.logger.setLevel(gunicorn_logger.level)

    root = logging.getLogger()
    if not os.path.exists('logs'):
        os.mkdir('logs')
    file_handler = RotatingFileHandler('logs/mine-server.log',
                                       maxBytes=1000000, backupCount=10)
    file_handler.setFormatter(logging.Formatter(
        '%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'))
    file_handler.setLevel(logging.DEBUG)
    root.addHandler(file_handler)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.WARNING)
    app.logger.addHandler(stream_handler)

    # Add any other (e.g. rdkit) logs to flask logs
    root.addHandler(default_handler)

    app.logger.setLevel(logging.DEBUG)
    app.logger.info('MINE-Server startup')
    app.logger.info('Running at http://127.0.0.1:5000')

    return app


if __name__ == "__main__":
    application = create_app()
    application.run(debug=False)
