# Import the core functionality and the app sub-packages
from . import core
from . import app

# Run the app locally at http://127.0.0.1:5000/
def launch_app(host="127.0.0.1", port=5000, debug=False):
    import webbrowser
    webbrowser.open("http://127.0.0.1:5000/")
    app.hes_off_app.run(host=host, port=port, debug=debug)