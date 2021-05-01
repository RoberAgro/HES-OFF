# Run the app locally at http://127.0.0.1:8080/
if __name__ == "__main__":
    from hes_off_app import app
    app.run(host="127.0.0.1", port=8080, debug=True)
