# Create app instance
from flask import Flask
hes_off_app = Flask(__name__)

# Import HES-OFF-app core functionality
from . import routes
from . import forms
from . import utilities


