
import os
import sys
import argparse
import uuid
import datetime
import requests
from lxml import etree
# import input_preprocessing as ip

# export SUMO_HOME="/usr/local/opt/sumo/share/sumo"

if 'SUMO_HOME' in os.environ:
    sumo_root = os.environ['SUMO_HOME']
    sumo_gui = sumo_root + "/bin/sumo-gui"
    sumo_emission = sumo_root + "/bin/emissionsDrivingCycle"
    sumo_phemlight = sumo_root + "/data/emissions/PHEMlight"
    tools = os.path.join(os.environ['SUMO_HOME'], 'tools')
    sys.path.append(tools)
else:
    sys.exit("please declare environment variable 'SUMO_HOME'")

import traci

def run():
    # step = 0
    while traci.simulation.getMinExpectedNumber() > 0:
        traci.simulationStep()
        # print(step)
        # step += 1
    traci.close()
    sys.stdout.flush()

def Main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', default="./emission-input/scenario2.sumocfg", type=str, help='path to sumo config file')
    parser.add_argument('--input', default="./emission-input/", type=str, help='path of emission inputs (= fcdoutput of sumo simulation)')
    parser.add_argument('--output', default="./emission-output/", type=str, help='path of emission outputs')
    args = parser.parse_args()
    sumo_emission_cmd = [sumo_emission, "-v", "-t", args.input + "trajoutput.xml", "--phemlight-path", sumo_phemlight, "--output", args.output + "emission_cycle_output.xml", "--emission-output", args.output + "emission_output.xml", "--sum-output", args.output + "emission_sum_output.xml"]
    cmd = ""
    for val in sumo_emission_cmd:
        cmd += val + " "
    print(cmd)
    os.system(cmd)
    # print(sumo_emission_cmd)
    # traci.start(sumo_emission_cmd)
    # run()

if __name__ == "__main__":
    print("running SUMO emissionsDrivingCycle tool")
    Main()