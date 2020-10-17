# Acinetobacter Dashboard
 [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0) [![Build Status](https://travis-ci.org/ba1/aci-dash.svg?branch=master)](https://travis-ci.org/ba1/aci-dash) [![Updates](https://pyup.io/repos/github/ba1/aci-dash/shield.svg)](https://pyup.io/repos/github/ba1/aci-dash/) [![Python 3](https://pyup.io/repos/github/ba1/aci-dash/python-3-shield.svg)](https://pyup.io/repos/github/ba1/aci-dash/) [![Coverage](https://codecov.io/github/ba1/aci-dash/coverage.svg?branch=master)](https://codecov.io/github/ba1/aci-dash?branch=master) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)


Integration of Comparative Genomics Analyses across the entire genus Acinetobacter.


## Usage
If you are on Linux and you have [virtualenvwrapper](https://virtualenvwrapper.readthedocs.io/en/latest/) installed you can run the script `utility/setup_virtualenv_and_repo.sh` to:

- create a python virtual environment and activate it
- install all project dependencies from `requirements.txt`
- create a git repository
- create your `Initial commit`

Here is how you run the script:

```shell
cd aci-dash
# mind the dot!
. utility/setup_virtualenv_and_repo.sh
```

Then you will need to create an `.env` file where to store your environment variables (SECRET key, plotly credentials, API keys, etc). Do NOT TRACK this `.env` file. See `.env.example`.

Run all tests with a simple:

```
pytest -v
```


## Run your Dash app
Check that the virtual environment is activated, then run:

```shell
cd aci_dash
python app.py
```

## Code formatting
To format all python files, run:

```shell
black .
```

## Pin your dependencies

```shell
pip freeze > requirements.txt
```

## Deploy on Heroku
Follow the [Dash deployment guide](https://dash.plot.ly/deployment) or have a look at the [dash-heroku-template](https://github.com/plotly/dash-heroku-template)
