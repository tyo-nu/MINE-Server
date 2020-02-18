import logging
import os
from logging.handlers import RotatingFileHandler

from api.config import Config
from flask import Flask
from flask.logging import default_handler
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

    # Initialize logger
    root = logging.getLogger()
    if not os.path.exists('logs'):
        os.mkdir('logs')
    file_handler = RotatingFileHandler('logs/mine-server.log',
                                       maxBytes=1000000, backupCount=10)
    file_handler.setFormatter(logging.Formatter(
        '%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'))
    file_handler.setLevel(logging.DEBUG)

    # Add other (e.g. cobrapy) logs to flask logs
    root.addHandler(default_handler)

    app.logger.setLevel(logging.DEBUG)
    app.logger.info('MINE-Server startup')

    if app.config['MODE'] == 'development':
        app.logger.info('Running in development mode at 127.0.0.1:5000')
    elif app.config['MODE'] == 'production':
        app.logger.info('Running in production mode.')
    elif app.config['MODE'] == 'testing':
        app.logger.info('Running in testing mode.')
    else:
        raise ValueError('MODE in config.py must be one of the following: ',
                         '"production", "development", or "testing".')

    return app


if __name__ == "__main__":
    app = create_app()
    app.run(debug=True)
