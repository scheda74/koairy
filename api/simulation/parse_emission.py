import os
import sys
import argparse
import uuid
import datetime
import requests
import xmltodict
import asyncio
import json
from lxml import etree

if 'SUMO_HOME' in os.environ:
    tools = os.path.join(os.environ['SUMO_HOME'], 'tools')
    sys.path.append(tools)
else:   
    sys.exit("please declare environment variable 'SUMO_HOME'")

import sumolib

net = None

async def parse_simulated_emissions():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', default=os.path.dirname(__file__) + "/data/emission-input/scenario2.sumocfg", type=str, help='path to sumo config file')
    parser.add_argument('--input', default=os.path.dirname(__file__) + "/data/emission-input/", type=str, help='path of emission inputs (= fcdoutput of sumo simulation)')
    parser.add_argument('--output', default=os.path.dirname(__file__) + "/data/emission-output/", type=str, help='path of emission outputs')
    args = parser.parse_args()
    
    global net
    net = sumolib.net.readNet(args.input + "kirchheim-street-names.net.xml")

    with open(args.output + "test_emission_output.xml") as fd:
        doc = xmltodict.parse(fd.read())
    
    emission_cycle = {}
    emission_sum = {}
    heatmap_data = {}
    max_co2 = 0
    for step in doc['emission-export']['timestep']:
        step_id = float(step["@time"])
        # print(step)
        # print("step_id: " + str(step_id))
        # emission_cycle[step_id] = {}
        emission_cycle[step_id] = {}
        emission_cycle[step_id]["step_id"] = step_id
        try:
            for vehicle in step["vehicle"]:
                veh_id = int(vehicle["@id"])
                veh_co2 = float(vehicle["@CO2"])
                veh_nox = vehicle["@NOx"] 
                veh_pmx = vehicle["@PMx"] # mg/s
                veh_fuel = vehicle["@fuel"] # ml/s
                lng, lat = net.convertXY2LonLat(float(vehicle["@x"]), float(vehicle["@y"]))

                coordinate = str(round(lat, 4)) + "," + str(round(lng, 4))
                # print(coordinate)
                if coordinate in heatmap_data:
                    heatmap_data[coordinate] += veh_co2
                else:
                    heatmap_data[coordinate] = veh_co2
                
                if heatmap_data[coordinate] > max_co2:
                    max_co2 = heatmap_data[coordinate]

                # print("veh_id: " + str(veh_id))
                if veh_id not in emission_cycle[step_id]:
                    emission_cycle[step_id][veh_id] = {}
                    emission_cycle[step_id][veh_id]["veh_id"] = veh_id
                
                emission_cycle[step_id][veh_id]["CO2"] = float(veh_co2)
                emission_cycle[step_id][veh_id]["NOx"] = float(veh_nox)
                emission_cycle[step_id][veh_id]["PMx"] = float(veh_pmx)
                emission_cycle[step_id][veh_id]["fuel"] = float(veh_fuel)
                emission_cycle[step_id][veh_id]["lat"] = lat
                emission_cycle[step_id][veh_id]["lng"] = lng

                if veh_id not in emission_sum:
                    emission_sum[veh_id] = {}
                    emission_sum[veh_id]["veh_id"] = veh_id
                    emission_sum[veh_id]["CO2"] = float(veh_co2)
                    emission_sum[veh_id]["NOx"] = float(veh_nox)
                    emission_sum[veh_id]["PMx"] = float(veh_pmx)
                    emission_sum[veh_id]["fuel"] = float(veh_fuel)
                else:
                    emission_sum[veh_id]["CO2"] += float(veh_co2)
                    emission_sum[veh_id]["NOx"] += float(veh_nox)
                    emission_sum[veh_id]["PMx"] += float(veh_pmx)
                    emission_sum[veh_id]["fuel"] += float(veh_fuel)
                
                # print(emission_sum[veh_id])
                # print(emission_cycle[step_id][veh_id])
            emission_cycle["veh_sum"] = emission_sum
        except (KeyError, TypeError) as err:
            print("Error: " + str(err))
            continue
    # print(heatmap_data)
    print("max CO2: " + str(max_co2))
    # with open(args.output + "emission_cycle.json", "w") as json_file:
    #     json.dump(heatmap_data, json_file)
    # return emission_cycle
    # WriteEmissionCycleToJson(emission_cycle, args.output + "emission_cycle.json")
    # print(emission_cycle)
    print(emission_sum)
    return heatmap_data

# if __name__ == "__main__":
#     print("running SUMO emissionsDrivingCycle tool")
#     Main()

