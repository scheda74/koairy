import os
import sys
import argparse
import uuid
import datetime
import requests
import xmltodict
import asyncio
import json
import xml.etree.ElementTree as ET
from simulation import calc_caqi as caqi
# import calc_caqi as caqi
# import xml.etree.cElementTree as etree
from lxml import etree


if 'SUMO_HOME' in os.environ:
    tools = os.path.join(os.environ['SUMO_HOME'], 'tools')
    sys.path.append(tools)
else:   
    sys.exit("please declare environment variable 'SUMO_HOME'")

import sumolib

net = None

# async def parse_simulated_emissions():
#     parser = argparse.ArgumentParser()
#     parser.add_argument('--config', default=os.path.dirname(__file__) + "/data/emission-input/scenario2.sumocfg", type=str, help='path to sumo config file')
#     parser.add_argument('--input', default=os.path.dirname(__file__) + "/data/emission-input/", type=str, help='path of emission inputs (= fcdoutput of sumo simulation)')
#     parser.add_argument('--output', default=os.path.dirname(__file__) + "/data/emission-output/", type=str, help='path of emission outputs')
#     args = parser.parse_args()
    
#     global net
#     net = sumolib.net.readNet(args.input + "kirchheim-street-names.net.xml")

#     with open(args.output + "emission_output.xml") as fd:
#         doc = xmltodict.parse(fd.read())
    
#     emission_cycle = {}
#     emission_sum = {}
#     heatmap_data = {}
#     max_co2 = 0
#     for step in doc['emission-export']['timestep']:
#         step_id = float(step["@time"])
#         # print(step)
#         # print("step_id: " + str(step_id))
#         # emission_cycle[step_id] = {}
#         emission_cycle[step_id] = {}
#         emission_cycle[step_id]["step_id"] = step_id
#         try:
#             for vehicle in step["vehicle"]:
#                 veh_id = int(vehicle["@id"])
#                 veh_co2 = float(vehicle["@CO2"])
#                 veh_nox = float(vehicle["@NOx"])
#                 veh_pmx = float(vehicle["@PMx"]) # mg/s
#                 veh_fuel = float(vehicle["@fuel"]) # ml/s
#                 lng, lat = net.convertXY2LonLat(float(vehicle["@x"]), float(vehicle["@y"]))

#                 coordinate = str(round(lat, 4)) + "," + str(round(lng, 4))
#                 # print(coordinate)
#                 if coordinate in heatmap_data:
#                     heatmap_data[coordinate] += veh_co2
#                 else:
#                     heatmap_data[coordinate] = veh_co2
                
#                 if heatmap_data[coordinate] > max_co2:
#                     max_co2 = heatmap_data[coordinate]

#                 # print("veh_id: " + str(veh_id))
#                 if veh_id not in emission_cycle[step_id]:
#                     emission_cycle[step_id][veh_id] = {}
#                     emission_cycle[step_id][veh_id]["veh_id"] = veh_id
                
#                 emission_cycle[step_id][veh_id]["CO2"] = float(veh_co2)
#                 emission_cycle[step_id][veh_id]["NOx"] = float(veh_nox)
#                 emission_cycle[step_id][veh_id]["PMx"] = float(veh_pmx)
#                 emission_cycle[step_id][veh_id]["fuel"] = float(veh_fuel)
#                 emission_cycle[step_id][veh_id]["lat"] = lat
#                 emission_cycle[step_id][veh_id]["lng"] = lng

#                 if veh_id not in emission_sum:
#                     emission_sum[veh_id] = {}
#                     emission_sum[veh_id]["veh_id"] = veh_id
#                     emission_sum[veh_id]["CO2"] = float(veh_co2)
#                     emission_sum[veh_id]["NOx"] = float(veh_nox)
#                     emission_sum[veh_id]["PMx"] = float(veh_pmx)
#                     emission_sum[veh_id]["fuel"] = float(veh_fuel)
#                 else:
#                     emission_sum[veh_id]["CO2"] += float(veh_co2)
#                     emission_sum[veh_id]["NOx"] += float(veh_nox)
#                     emission_sum[veh_id]["PMx"] += float(veh_pmx)
#                     emission_sum[veh_id]["fuel"] += float(veh_fuel)
                
#                 # print(emission_sum[veh_id])
#                 # print(emission_cycle[step_id][veh_id])
#             emission_cycle["veh_sum"] = emission_sum
#         except (KeyError, TypeError) as err:
#             print("Error: " + str(err))
#             continue
#     # print(heatmap_data)
#     print("max CO2: " + str(max_co2))
#     # with open(args.output + "emission_cycle.json", "w") as json_file:
#     #     json.dump(heatmap_data, json_file)
#     # return emission_cycle
#     # WriteEmissionCycleToJson(emission_cycle, args.output + "emission_cycle.json")
#     # print(emission_cycle)
#     print(emission_sum)
#     return heatmap_data

async def parse_emissions():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', default=os.path.dirname(__file__) + "/data/emission-input/scenario2.sumocfg", type=str, help='path to sumo config file')
    parser.add_argument('--input', default=os.path.dirname(__file__) + "/data/emission-input/", type=str, help='path of emission inputs (= fcdoutput of sumo simulation)')
    parser.add_argument('--output', default=os.path.dirname(__file__) + "/data/emission-output/", type=str, help='path of emission outputs')
    args = parser.parse_args()

    # net needed to convert x, z to lat, long 
    global net
    net = sumolib.net.readNet(args.input + "kirchheim-street-names.net.xml")

    # root = etree.parse(args.output + "test_emission_output.xml")
    xml_file = args.output + "emission_output.xml"
    emission_cycle = {}
    emission_sum = {}
    raw_pollution_data = {}
    max_co2 = 0
    # root = etree.Element("root") 
    for _, timestep in etree.iterparse(xml_file, tag="timestep"):
        step_id = float(timestep.attrib["time"])
        print("* step - %s -- %s" % (step_id, (step_id / 11800) * 100))
        # emission_cycle[step_id] = {}
        # emission_cycle[step_id]["step_id"] = step_id
        for vehicle in timestep.iterfind("vehicle"):
            veh_id = int(vehicle.attrib["id"])
            veh_noise = float(vehicle.attrib["noise"])
            veh_type = vehicle.attrib["type"]
            veh_eclass = vehicle.attrib["eclass"]
            veh_co2 = float(vehicle.attrib["CO2"])
            veh_co = float(vehicle.attrib["CO"])
            veh_hc = float(vehicle.attrib["HC"])
            veh_speed = float(vehicle.attrib["speed"])
            veh_nox = float(vehicle.attrib["NOx"]) 
            veh_pmx = float(vehicle.attrib["PMx"]) # mg/s
            veh_fuel = float(vehicle.attrib["fuel"]) # ml/s
            lng, lat = net.convertXY2LonLat(float(vehicle.attrib["x"]), float(vehicle.attrib["y"]))
            # print(veh_id)
            # print("** veh - %s" % veh_id)
            # print("*** - calculating coordinates")

            coordinate = str(round(lat, 4)) + "," + str(round(lng, 4))
            if coordinate in raw_pollution_data:
                raw_pollution_data[coordinate]["CO2"] += veh_co2
                raw_pollution_data[coordinate]["NOX"] += veh_nox
                raw_pollution_data[coordinate]["PMX"] += veh_pmx
                raw_pollution_data[coordinate]["CO"] = veh_co
            else:
                raw_pollution_data[coordinate] = {}
                raw_pollution_data[coordinate]["CO2"] = veh_co2
                raw_pollution_data[coordinate]["NOX"] = veh_nox
                raw_pollution_data[coordinate]["PMX"] = veh_pmx
                raw_pollution_data[coordinate]["CO"] = veh_co
            
            if raw_pollution_data[coordinate]["CO2"] > max_co2:
                max_co2 = raw_pollution_data[coordinate]["CO2"]
            
            # results = root.xpath("//location[@coord = '%s']" % coordinate)
            # if results:
            #     location = results[0]
            #     location.attrib["value"] += str(veh_co2)
            # else:
            #     child = etree.SubElement(root, "location")
            #     child.attrib["coord"] = coordinate
            #     child.attrib["value"] = str(veh_co2)
            

            # print("veh_id: " + str(veh_id))
            # if veh_id not in emission_cycle[step_id]:
            #     emission_cycle[step_id][veh_id] = {}
            #     emission_cycle[step_id][veh_id]["veh_id"] = veh_id
            
            # emission_cycle[step_id][veh_id]["CO2"] = float(veh_co2)
            # emission_cycle[step_id][veh_id]["NOx"] = float(veh_nox)
            # emission_cycle[step_id][veh_id]["PMx"] = float(veh_pmx)
            # emission_cycle[step_id][veh_id]["fuel"] = float(veh_fuel)
            # emission_cycle[step_id][veh_id]["lat"] = lat
            # emission_cycle[step_id][veh_id]["lng"] = lng

            # if veh_id not in emission_sum:
            #     emission_sum[veh_id] = {}
            #     emission_sum[veh_id]["veh_id"] = veh_id
            #     emission_sum[veh_id]["CO2"] = float(veh_co2)
            #     emission_sum[veh_id]["NOx"] = float(veh_nox)
            #     emission_sum[veh_id]["PMx"] = float(veh_pmx)
            #     emission_sum[veh_id]["fuel"] = float(veh_fuel)
            # else:
            #     emission_sum[veh_id]["CO2"] += float(veh_co2)
            #     emission_sum[veh_id]["NOx"] += float(veh_nox)
            #     emission_sum[veh_id]["PMx"] += float(veh_pmx)
            #     emission_sum[veh_id]["fuel"] += float(veh_fuel)
            # print(etree.tostring(vehicle))
            # print("\n\n\n\n\n\n")
        
        # emission_cycle["veh_sum"] = emission_sum
        # print('%s' % (timestep.findtext('vehicle')))
        # print(etree.tostring(timestep))
        # print("\n\n\n\n\n\n")
        timestep.clear(keep_tail=True)
    print(raw_pollution_data)
    return raw_pollution_data
    # return caqi.calc_sub_indices(raw_pollution_data)
    # return calc_sub_indices(raw_pollution_data)
    # with open(args.output + "emission_cycle.json", "w") as json_file:
    #     json.dump(heatmap_data, json_file)
    # return "yes"
    # return heatmap_data

async def get_caqi_data():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default=os.path.dirname(__file__) + "/data/emission-output/", type=str, help='path of emission outputs')
    args = parser.parse_args()

    raw_pollution_data = await parse_emissions()
    air = caqi.calc_indices(raw_pollution_data)
    with open(args.output + "temp_caqi.json", "w") as json_file:
        json.dump(air, json_file)
    return air

