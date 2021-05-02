# HES-OFF
## Description

A research-based prototype tool in Python

## Installation

It is planned to add the `hes_off` package to the `pip` and `conda` package managers once the code is more mature. For the time being, it is possible to import the local `hes-off` package if your python script is within the HES-OFF project.

If you want to import the `hes_off` package from a python script located anywhere in you computer, you can add the location of the `HES-OFF` project to the `PYTHONPATH` environmental variable. Here are [some instructions](https://stackoverflow.com/questions/3701646/how-to-add-to-the-pythonpath-in-windows-so-it-finds-my-modules-packages) about how to create and modify environmental variables in Windows OS. In my case, the project was located at `C:\Users\rober\git\HES-OFF` 

## Getting started



## Documentation

It is planned to create a comprehensive documentation page using [Read the Docs.](https://readthedocs.org/)



## Related publications



## License

The `hes_off` package is licensed under the terms of the MIT license. See the [license file](LICENSE.md) for more information.





# Requirements
Create a Conda environment from .yml to install dependencies

​	conda env create environment.yml



Create a list of requirements for the Heroku deployment

​		pip list --format=freeze > requirements.txt





Create Heroku account online

https://dashboard.heroku.com/apps



Install the command line interface (CLI):

https://devcenter.heroku.com/articles/heroku-cli



Create "Procfile" in the app root directory:

```shell
echo "web: gunicorn app:app" > Procfile
```



Push the heroku app

```
git push heroku main
```







https://stackoverflow.com/questions/47304291/heroku-upload-could-not-find-a-version-that-satisfies-the-requirement-anaconda/56754565

Remove items in requirements.txt if necessary




