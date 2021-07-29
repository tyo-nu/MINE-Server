""" Module for all Exception subclassses."""


class InvalidUsage(Exception):
    """Custom Exception handler for the API. See source at
    http://flask.pocoo.org/docs/1.0/patterns/apierrors/ for more info.

    Attributes
    ----------
    message : str
        Human readable string describing the exception.
    status_code : int, optional
        HTTP status code - defaults to 400 (Bad Request).
    payload : dict, optional
        Extra information on the error.
    """

    def __init__(self, message, status_code=400, payload=None):
        Exception.__init__(self)
        self.message = message
        if status_code is not None:
            self.status_code = status_code
        self.payload = payload

    def to_dict(self):
        """Convert payload to dict with message."""
        rv = dict(self.payload or ())
        rv['message'] = self.message
        return rv
