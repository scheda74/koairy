import os, sys

if 'SUMO_HOME' in os.environ:
    tools = os.path.join(os.environ['SUMO_HOME'], 'tools')
    sys.path.append(tools)
else:
    sys.exit("please declare environment variable 'SUMO_HOME'")

class Constants:
    BASEDIR = os.path.dirname(__file__)
    DEFAULT_NET_INPUT = BASEDIR + "/data/traffic-input/road-network-default.net.xml"
    AREA_OF_INTEREST = BASEDIR + "/data/traffic-input/areas-of-interest.taz.xml"
    TRIP_OUTPUT = BASEDIR + "/data/traffic-input/trip-"
    ROUTE_OUTPUT = BASEDIR + "/data/traffic-input/route-"
    WEIGHT_INPUT = BASEDIR + "/data/traffic-input/weight-"
    EMISSION_OUTPUT_BASE = BASEDIR + "/data/emission-output/"
    EMISSION_OUTPUT = BASEDIR + "/data/emission-output/emission_output.xml"
    SUMO_CFG = BASEDIR + "/data/traffic-input/simulation-"

    SUMO_ROOT = os.environ['SUMO_HOME']
    SUMO_GUI = SUMO_ROOT + "/bin/sumo-gui"
    SUMO_COMMANDLINE = SUMO_ROOT + "/bin/sumo"
    SUMO_EM_DRIVING_CYCLE = SUMO_ROOT + "/bin/emissionsDrivingCycle"
    RANDOM_TRIP_TOOL = tools + "/randomTrips.py"

    PHEMLIGHT_PATH = SUMO_ROOT + "/data/emissions/PHEMlight"
    
    VALID_AREA_IDS = {
        'aschheim_west',
        'ebersberg_east',
        'feldkirchen_west',
        'heimstetten_industrial_1',
        'heimstetten_industrial_2',
        'heimstetten_residential',
        'kirchheim_industrial_east',
        'kirchheim_industrial_west',
        'kirchheim_residential',
        'unassigned_edges'
    }