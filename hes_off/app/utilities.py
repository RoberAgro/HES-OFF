import io
import base64
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.backends.backend_svg import FigureCanvasSVG


def file_allowed(filename, allowed_extensions):
    """Does filename have the right extension?"""
    return "." in filename and filename.rsplit('.', 1)[1] in allowed_extensions

def convert_fig_to_png(fig):
    """Convert Matplotlib fig to .png string"""
    pngImage = io.BytesIO()
    FigureCanvasAgg(fig).print_png(pngImage)
    pngImageString = "data:image/png;base64,"
    pngImageString += base64.b64encode(pngImage.getvalue()).decode('utf8')
    return pngImageString

def convert_fig_to_svg(fig):
    """Convert Matplotlib fig to .svg string"""
    svgImage = io.BytesIO()
    FigureCanvasSVG(fig).print_svg(svgImage)
    svgImageString = "data:image/svg+xml;base64,"
    svgImageString += base64.b64encode(svgImage.getvalue()).decode('utf8')
    return svgImageString

def update_form_from_dict(form, field_dict):
    """Update WTForm field values from Python dictionary"""
    for field in form:
        if field.name in field_dict.keys():
            value = field_dict[field.name]
            if field.type == "FieldList":
                if not isinstance(value, list):
                    value = [value]
                for index in range(len(field)):
                    try:
                        field[index].data = value[index]
                        field[index].raw_data = None
                    except:
                        field[index].data = None
                        field[index].raw_data = None
            else:
                field.data = value
                field.raw_data = None

    return form


def create_dict_from_form(form):
    """Create Python dictionary from WTForm field values"""

    # Instantiate dictionary
    field_dict = {}

    # Fill dictionary with form entries
    for fieldname, value in form.data.items():
        field_dict[fieldname] = value

    # Convert "decimal" entries to "float"
    for key, value in field_dict.items():
        if isinstance(value, list):
            try:
                value = [float(item) for item in value]
            except:
                pass
        else:
            try:
                value = float(value)
            except:
                pass

        field_dict[key] = value

    return field_dict


def initialize_plot_dict():
    """Initialize a Python dictionary containing the plots for the UI"""
    plot_dict = {}
    plot_dict["CO2_EMISSIONS"] = ["CO<sub>2</sub>&nbsp;emissions graph", ""]
    plot_dict["H2_LEVEL"]      = ["Hydrogen level graph", ""]
    plot_dict["POWER_DEFICIT"] = ["Power deficit graph", ""]
    # plot_dict["WIND_SPEED"]    = ["Wind speed graph", ""]
    # plot_dict["POWER_BALANCE"]    = ["Power balance graph", ""]
    return plot_dict


def create_plot_dict(EnergySystem):
    """Create a Python dictionary containing the plots for the UI"""
    plot_dict = initialize_plot_dict()
    plot_dict["CO2_EMISSIONS"][-1] = EnergySystem.plot_carbon_dioxide_emissions()
    plot_dict["H2_LEVEL"][-1]      = EnergySystem.plot_hydrogen_level()
    plot_dict["POWER_DEFICIT"][-1] = EnergySystem.plot_power_deficit()
    # plot_dict["WIND_SPEED"][-1]    = EnergySystem.plot_wind_timeseries()
    # plot_dict["POWER_BALANCE"][-1]    = EnergySystem.plot_power_balance()
    plot_dict = {key:[label,convert_fig_to_svg(value)] for key, (label, value) in plot_dict.items()}
    return plot_dict


def initialize_data_dict():
    """Initialize a Python dictionary containing the data for the UI"""
    data_dict = {}
    data_dict["HEADER"] = ["Performance indicator", "", ["Peak years", "Mid years", "Tail years", "Total"]]
    data_dict["H2_utilized"] = ["H<sub>2</sub>&nbsp;consumption (ton)", "Accumulated hydrogen consumption", 4 * [""]]
    # data_dict["NG_utilized"] = ["NG consumption (kton)", "Accumulated natural gas consumption", 4 * [""]]
    data_dict["CO2_emissions"] = ["CO<sub>2</sub>&nbsp;emissions (kton)", "Accumulated carbon dioxide emissions", 4*[""]]
    data_dict["GT_energy"] = ["GT energy  (GWe&#183h)", "Accumulated gas turbine electricity generation", 4 * [""]]
    data_dict["WT_energy"] = ["WT energy (GWe&#183h)", "Accumulated wind farm electricity generation", 4 * [""]]
    data_dict["FC_energy"] = ["FC energy (GWe&#183h)", "Accumulated fuel cell electricity generation", 4 * [""]]
    data_dict["EL_energy"] = ["EL energy (GWe&#183h)", "Accumulated electrolyzer electricity consumption", 4 * [""]]
    data_dict["energy_deficit"] = ["Energy deficit (GWe&#183h)", "Accumulated energy deficit", 4*[""]]
    return data_dict


def create_data_dict(EnergySystem):
    """Create a Python dictionary containing the data for the UI"""
    data_dict = initialize_data_dict()
    data_dict["H2_utilized"][-1]    = ["{:0.2f}".format(item*3600/1e3) for item in EnergySystem.H2_utilized] # Convert kg/s to ton
    # data_dict["NG_utilized"][-1]    = ["{:0.2f}".format(item*3600/1e6) for item in EnergySystem.NG_utilized] # Convert kg/s to kton
    data_dict["CO2_emissions"][-1]  = ["{:0.2f}".format(item/1e6) for item in EnergySystem.CO2_emissions] # Convert kg to kton
    data_dict["GT_energy"][-1]      = ["{:0.2f}".format(item/1e9) for item in EnergySystem.GT_energy]
    data_dict["WT_energy"][-1]      = ["{:0.2f}".format(item/1e9) for item in EnergySystem.WT_energy]
    data_dict["FC_energy"][-1]      = ["{:0.2f}".format(item/1e9) for item in EnergySystem.FC_energy]
    data_dict["EL_energy"][-1]      = ["{:0.2f}".format(item/1e9) for item in EnergySystem.EL_energy]
    data_dict["energy_deficit"][-1] = ["{:0.2f}".format(max(0.0, item/1e9)) for item in EnergySystem.energy_deficit]
    return data_dict
