# Run the app locally at http://127.0.0.1:8080/
from hes_off_app import app
if __name__ == "__main__":
    app.run(host="127.0.0.1", port=8080, debug=True)
