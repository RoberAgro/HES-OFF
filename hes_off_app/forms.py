# Import packages
import flask_wtf
import wtforms as wtf


# Define placeholders
placeholder_positive = "Positive value"
placeholder_integer = "Positive integer"
placeholder_percentage = "Percentage"


# Define custom validators
def IntegerRequired():
    message = "Must be an integer value"

    def _IntegerRequired(form, field):
        try:
            float(field.data)
        except:
            raise wtf.ValidationError(message)

        if not float(field.data).is_integer():
            raise wtf.ValidationError(message)

    return _IntegerRequired

def RechargeCofireTest():
    message = "Hydrogen recharge threshold must be lower than the co-fire threshold"

    def _RechargeCofireTest(form, field):

        try:
            H2_COFIRE_THRESHOLD = float(form.H2_COFIRE_THRESHOLD.data)
            H2_RECHARGE_THRESHOLD = float(form.H2_RECHARGE_THRESHOLD.data)
        except:
            raise wtf.ValidationError("Something went wrong")

        if H2_RECHARGE_THRESHOLD > H2_COFIRE_THRESHOLD:
            raise wtf.ValidationError(message)

    return _RechargeCofireTest


# Define the form
class HesOffForm(flask_wtf.FlaskForm):

    # Disable CSRF
    class Meta:
        csrf = False

    # UI widgets
    INPUT_FILE = wtf.FileField(
        label="Load input from file",
        validators=[],
        render_kw={"data-bs-toggle": "tooltip",
                   "data-bs-placement": "top",
                   "title": "Select a configuration file containing the model parameters"})

    UPLOAD = wtf.SubmitField(
        label="Upload file",
        validators=[],
        render_kw={"formnovalidate": "formnovalidate",
                   "data-bs-toggle": "tooltip",
                   "data-bs-placement": "top",
                   "title": "Upload a configuration file with the model parameters"})

    DEFAULT = wtf.SubmitField(
        label="Default inputs",
        render_kw={"formnovalidate": "formnovalidate",
                   "data-bs-toggle": "tooltip",
                   "data-bs-placement": "top",
                   "title": "Set default model parameters"})

    COMPUTE = wtf.SubmitField(
        label="Compute results",
        render_kw={"data-bs-toggle": "tooltip",
                   "data-bs-placement": "top",
                   "title": "Evaluate the model and display the results"})

    PLOTS = wtf.SubmitField(
        label="Plots results",
        render_kw={"data-bs-toggle": "tooltip",
                   "data-bs-placement": "top",
                   "title": "Evaluate the model and visualize the results"})

    # Process specifications
    HEAT_DEMAND = wtf.FieldList(wtf.DecimalField(
        label="Heat demand (MW)",
        render_kw={"data-bs-toggle": "tooltip",
                   "data-bs-placement": "right",
                   "title": "Specify the process heat demand over the life-time of the platform"},
        default=None, places=2,
        validators=[wtf.validators.InputRequired(), wtf.validators.NumberRange(min=0)]),
        min_entries=3, max_entries=3)

    POWER_DEMAND = wtf.FieldList(wtf.DecimalField(
        label="Power demand (MW)",
        render_kw={"data-bs-toggle": "tooltip",
                   "data-bs-placement": "right",
                   "title": "Specify the process power demand over the life-time of the platform"},
        default=None, places=2,
        validators=[wtf.validators.InputRequired(), wtf.validators.NumberRange(min=0)]),
        min_entries=3, max_entries=3)

    STAGE_LENGTH = wtf.FieldList(wtf.DecimalField(
        label="Stage length (years)",
        render_kw={"data-bs-toggle": "tooltip",
                   "data-bs-placement": "right",
                   "title": "Specify the length of the peak, midlife and tail stages of the platform"},
        default=None, places=2,
        validators=[wtf.validators.InputRequired(), wtf.validators.NumberRange(min=0)]),
        min_entries=3, max_entries=3)

    # Gas turbine specifications
    GT_MODEL = wtf.SelectField(
        label="Gas turbine model",
        render_kw={"data-bs-toggle": "tooltip",
                   "data-bs-placement": "right",
                   "title": "Specify the model of the gas turbines"},
        choices=[("", ""), ("LM2500+G4", "LM2500+G4 (32 MW)"), ("LM6000-PF", "LM6000-PF (42 MW)")],
        validators=[wtf.validators.InputRequired()])

    GT_UNITS = wtf.DecimalField(
        label="Gas turbine units",
        render_kw={"placeholder": placeholder_integer,
                   "data-bs-toggle": "tooltip",
                   "data-bs-placement": "right",
                   "title": "Specify the number of gas turbines units"},
        default=None, places=0,
        validators=[wtf.validators.InputRequired(), wtf.validators.NumberRange(min=1), IntegerRequired()])

    GT_MAX_H2 = wtf.DecimalField(
        label="Max. hydrogen content (%)",
        render_kw={"placeholder": placeholder_percentage,
                   "data-bs-toggle": "tooltip",
                   "data-bs-placement": "right",
                   "title": "Specify the maximum volumetric fraction of hydrogen allowed in the gas turbine fuel mixture"},
        default=None, places=2,
        validators=[wtf.validators.InputRequired(), wtf.validators.NumberRange(min=0, max=100)])

    # Wind turbine specifications
    WT_MODEL = wtf.SelectField(
        label="Wind turbine model",
        render_kw={"data-bs-toggle": "tooltip",
                   "data-bs-placement": "right",
                   "title": "Specify the model of the wind turbines"},
        choices=[("", ""), ("HYWIND", "Hywind"), ("NREL", "NREL")],
        validators=[wtf.validators.InputRequired()])

    WT_RATED_POWER = wtf.DecimalField(
        label="Wind farm rated power (MW)",
        render_kw={"placeholder": placeholder_positive,
                   "data-bs-toggle": "tooltip",
                   "data-bs-placement": "right",
                   "title": "Specify the rated power of the wind farm"},
        default = None, places = 2,
        validators=[wtf.validators.InputRequired(), wtf.validators.NumberRange(min=0)])

    WT_HUB_HEIGHT = wtf.DecimalField(
        label="Wind turbine hub height (m)",
        render_kw={"placeholder": placeholder_positive,
                   "data-bs-toggle": "tooltip",
                   "data-bs-placement": "right",
                   "title": "Specify the height of the wind turbine hub (used to compute the wind speed factor)"},
        default=None, places=2,
        validators=[wtf.validators.InputRequired(), wtf.validators.NumberRange(min=0)])

    WT_REF_HEIGHT = wtf.DecimalField(
        label="Wind data reference height (m)",
        render_kw={"placeholder": placeholder_positive,
                   "data-bs-toggle": "tooltip",
                   "data-bs-placement": "right",
                   "title": "Specify the height of the wind speed data (used to compute the wind speed factor)"},
        default=None, places=2,
        validators=[wtf.validators.InputRequired(), wtf.validators.NumberRange(min=0)])

    WIND_FILENAME = wtf.SelectField(
        label="Wind data file",
        render_kw={"data-bs-toggle": "tooltip",
                   "data-bs-placement": "right",
                   "title": "Specify the name of the file containing the wind speed data"},
        choices=[("", ""), ("SLEIPNERWIND", "Sleipner wind")],
        validators=[wtf.validators.InputRequired()])

    # Electrolyzer specifications
    EL_MODEL = wtf.SelectField(
        label="Electrolyzer model",
        render_kw={"data-bs-toggle": "tooltip",
                   "data-bs-placement": "right",
                   "title": "Specify the model of the electrolyzer system"},
        choices=[("", ""), ("NEL_HYDROGEN", "Nel Hydrogen")],
        # choices=["", "NEL_HYDROGEN", "POLYNOMIAL_EFFICIENCY"],
        validators=[wtf.validators.InputRequired()])

    EL_RATED_POWER = wtf.DecimalField(
        label="Electrolyzer rated power (MW)",
        render_kw={"placeholder": placeholder_positive,
                   "data-bs-toggle": "tooltip",
                   "data-bs-placement": "right",
                   "title": "Specify the rated power of the electrolyzer system"},
        default=None, places=2,
        validators=[wtf.validators.InputRequired(), wtf.validators.NumberRange(min=0)])

    # EL_EFFICIENCY = wtf.FieldList(wtf.DecimalField(
    #     label="Electrolyzer efficiency (polynomial coefficients)",
    #     render_kw={"placeholder": "Scalar value", "disabled": ""},
    #     default=None, places=4,
    #     validators=[]), min_entries=4, max_entries=5)

    # Fuel cell specifications
    FC_MODEL = wtf.SelectField(
        label="Fuel cell model",
        render_kw={"data-bs-toggle": "tooltip",
                   "data-bs-placement": "right",
                   "title": "Specify the model of the fuel cell system"},
        choices=[("", ""), ("POWERCELL_S3", "PowerCell S3")],
        # choices=["", "POWERCELL_S3", "POLYNOMIAL_EFFICIENCY"],
        validators=[wtf.validators.InputRequired()])

    FC_RATED_POWER = wtf.DecimalField(
        label="Fuel cell rated power (MW)",
        render_kw={"placeholder": placeholder_positive,
                   "data-bs-toggle": "tooltip",
                   "data-bs-placement": "right",
                   "title": "Specify the rated power of the fuel cell system"},
        default=None, places=2,
        validators=[wtf.validators.InputRequired(), wtf.validators.NumberRange(min=0)])

    # FC_EFFICIENCY = wtf.FieldList(wtf.DecimalField(
    #     label="Fuel cell efficiency (polynomial coefficients)",
    #     default=None, places=4,
    #     render_kw={"placeholder": "Scalar value", "disabled": ""},
    #     validators=[]), min_entries=4, max_entries=5)

    # Hydrogen storage specifications
    H2_CAPACITY = wtf.DecimalField(
        label="Hydrogen storage capacity (kg)",
        render_kw={"placeholder": placeholder_positive,
                   "data-bs-toggle": "tooltip",
                   "data-bs-placement": "right",
                   "title": "Specify the total capacity of the hydrogen storage system"},
        default=None, places=1,
        validators=[wtf.validators.InputRequired(), wtf.validators.NumberRange(min=0)])

    H2_INITIAL_LEVEL = wtf.DecimalField(
        label="Initial level of hydrogen (%)",
        render_kw={"placeholder": placeholder_percentage,
                   "data-bs-toggle": "tooltip",
                   "data-bs-placement": "right",
                   "title": "Specify the initial level of the hydrogen storage system"},
        default=None, places=2,
        validators=[wtf.validators.InputRequired(), wtf.validators.NumberRange(min=0, max=100)])

    H2_RECHARGE_THRESHOLD = wtf.DecimalField(
        label="Hydrogen recharge threshold (%)",
        render_kw={"placeholder": placeholder_percentage,
                   "data-bs-toggle": "tooltip",
                   "data-bs-placement": "right",
                   "title": "Specify the storage level below which the gas turbine is used to produce hydrogen"},
        default=None, places=2,
        validators=[wtf.validators.InputRequired(), wtf.validators.NumberRange(min=0, max=100)])

    H2_COFIRE_THRESHOLD = wtf.DecimalField(
        label="Hydrogen co-fire threshold (%)",
        render_kw={"placeholder": placeholder_percentage,
                   "data-bs-toggle": "tooltip",
                   "data-bs-placement": "right",
                   "title": "Specify the storage level above which hydrogen is co-fired in the gas turbines"},
        default=None, places=2,
        validators=[wtf.validators.InputRequired(), wtf.validators.NumberRange(min=0, max=100), RechargeCofireTest()])
