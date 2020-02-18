"""Used just to create the mongo connection. See Stack Overflow question
33166612 for details on why this needs to be its own file."""

from flask_pymongo import PyMongo


mongo = PyMongo()
