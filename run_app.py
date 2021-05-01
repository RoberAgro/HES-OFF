# Load the app
from hes_off_app import app

# Run the app locally at http://127.0.0.1:5000/
if __name__ == "__main__":
    app.run(host="127.0.0.1", port=5000, debug=True)
