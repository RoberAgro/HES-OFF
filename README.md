# HES-OFF
## Description

Despite the efforts to move towards a carbon-neutral future, the global demand of fossil fuels continues to grow.  As a result, the oil and gas (O&G) industry if facing a double challenge: meet the increased need for energy while simultaneously reducing its overall emissions. One way to achieve these goals is through the use of hybrid energy systems that combine renewable energies with conventional power generation.

In this context, `HES-OFF` Python package offers a platform-independent environment for the simulation of hybrid energy systems for stable heat and power supply in off-shore oil and gas installations. The code was developed aiming at ease of use and computational speed. Some notable features of the package are listed below:

- Supports the simulation of off-shore energy systems that integrate offshore wind power with gas turbines and an energy storage solution based on proton exchange membrane fuel cells and electrolysers.

- Contains commercial and user-defined numerical models to simulate the performance of the different system components.

- Includes an intuitive, web-based graphical user interface to perform the simulations and visualize the results.
  - The backend of the interface was developed using the Flask framework
  - The frontend of the interface was developed using the Bootstrap framework
  - The web-based app is currently hosted in the Heroku platform
  
- Uses the Numba JIT compiler to achieve execution times comparable with those of compiled languages such as C++ or FORTRAN.

  

## Getting started

#### Running the all online

There is not need to install anything to start using the online web application! Just [open the app](https://hes-off.herokuapp.com/) and start playing around with the different parameters and settings. The user interface has many tooltips and warnings to guide the users and help them understand what inputs are expected.

#### Running the app locally

If you want to have access to the core functionality of the HES-OFF package, you can install it with [pip](https://pip.pypa.io/en/stable/): 

```bash
pip install hes_off
```

Once the package is installed, you can run the web application locally typing the following command in your terminal:

```shell
python -c "import hes_off; hes_off.launch_app()"
```



## Documentation

To be completed

- [ ] Add link to Read the Docs documentation
- [ ] Add the links to related publications



## License

The `hes_off` package is licensed under the terms of the MIT license. See the [license file](LICENSE.md) for more information.



## Contact information

To be completed



## To-do list

- [ ] Add option to export Excel files with results (Roberto)
- [ ] Add option to add a user-defined model for the electrolyzer and fuel cell systems (Roberto)
- [ ] Prepare some test cases and compare results with MATLAB (Luca)
- [ ] Prepare a documentation page using  [Read the Docs](https://readthedocs.org/) (Erick)
- [ ] Implement the electric grid surrogate model (Erick)
- [ ] Implement a gradient-free optimization for design
- [ ] Prepare a paper for an open-software journal
  - [Journal of open-source software](https://joss.theoj.org/)
  - [Elsevier's SoftwareX journal](https://www.journals.elsevier.com/softwarex/)



## Pip deployment

_This section is intended for developers_

In order to add the `hes_off` package to the pip installer, you must create a PyPi account. Once you are set, open a terminal in the project's root directory and type the following commands:

1. Remove any existing `build`, `dist` and `egg` folders:

   ```shell
   rm -r build dist *.egg-info
   ```

2. Create a source archive and a wheel for your package:

    ```shell
    python3 setup.py sdist bdist_wheel
    ```

3. Upload the package to PyPi using twine (you may need to install twine):

    ```shell
    python -m twine upload dist/*
    ```
    
4. Fill in your user name and credentials and wait until the package upload is complete



## Heroku deployment

*This section is intended for developers*

In order to deploy the `hes_off` app, you must create a [Heroku account](https://dashboard.heroku.com/apps) and install the [Heroku command line interface](https://devcenter.heroku.com/articles/heroku-cli). Once you are set, open a terminal in the project's root directory and type the following commands:

1. Create a Conda virtual environment with all the dependencies using the `environment.yml` file:

   ```shell
   conda env create environment.yml
   ```

2. Check that the app works locally typing this one-liner in your terminal:

   ```shell
   python -c "import hes_off; hes_off.launch_app()"
   ```

3. If the app works locally, then create a `requirements.txt` file with the list of dependencies:

   ```shell
   pip list --format=freeze > requirements.txt
   ```

4. Create the `Procfile` file in case it does not exists already:

   ```shell
   echo "web: gunicorn app:app" > Procfile
   ```

5. Deploy the app to Heroku with a Git push;

   ```shell
   git push heroku main
   ```

That's it! The app might take some moments to go online.

**Known caveats**

In some cases the `requirements.txt` file might include some unnecessary packages that cause the Heroku deployment to fail. If this is the case, [remove the problematic packages](https://stackoverflow.com/questions/47304291/heroku-upload-could-not-find-a-version-that-satisfies-the-requirement-anaconda/56754565) and try to re-deploy.


