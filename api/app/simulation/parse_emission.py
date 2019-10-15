
import os
import sys
import json
import asyncio
import pandas as pd
import numpy as np

from lxml import etree
from operator import itemgetter
from timeit import default_timer as timer
from app.simulation import calc_caqi as caqi
from app.simulation import preprocessor as ip
from app.simulation.constants import Constants

if 'SUMO_HOME' in os.environ:
    tools = os.path.join(os.environ['SUMO_HOME'], 'tools')
    sys.path.append(tools)
else:   
    sys.exit("please declare environment variable 'SUMO_HOME'")

import sumolib
net = sumolib.net.readNet(Constants.DEFAULT_NET_INPUT)

def extract_attributes(context, fields):
    values = itemgetter(*fields)
    for _, elem in context:
        yield values(elem.attrib)
        elem.clear()
        while elem.getprevious() is not None:
            del elem.getparent()[0]
    del context

def parse_emissions(filename):
    context = etree.iterparse(filename, tag="vehicle")

    # create a dataframe from XML data a single call
    coords = ['x', 'y']
    entries = ['CO2', 'CO', 'NOx', 'PMx', 'fuel']
    df = pd.DataFrame(
        extract_attributes(context, coords + entries),
        columns=coords + entries, dtype=np.float)

    # convert *all coordinates together*, remove the x, y columns
    # note that the net.convertXY2LonLat() call *alters the 
    # numpy arrays in-place* so we don’t want to keep them anyway. 
    df['lng'], df['lat'] = net.convertXY2LonLat(df.x.to_numpy(), df.y.to_numpy())
    df.drop(coords, axis=1, inplace=True)

    # 'group' data by rounding the latitude and longitude
    # effectively creating areas of 1/10000th degrees per side
    latlng = ['lat', 'lng']
    df[latlng] = df[latlng].round(4)

    # aggregate the results and return summed dataframe
    return df.groupby(latlng)[entries].sum()

def get_caqi_data():
    timer_start = timer()

    emissions = parse_emissions(Constants.EMISSION_OUTPUT)
    seconds = timer() - timer_start
    print(emissions)
    print("[etree] Finished parsing XML in %s seconds" % seconds)
    caqi_emissions = emissions.apply(caqi.calc_indices, axis=1)
    # print(caqi_emissions)
    return caqi_emissions.to_json(orient='index')
    # df = pd.DataFrame(emissions)
    # df.app
    # return emissions



# def init_env():
#     parser = argparse.ArgumentParser()
#     parser.add_argument('--config', default=os.path.dirname(__file__) + "/data/emission-input/scenario2.sumocfg", type=str, help='path to sumo config file')
#     parser.add_argument('--input', default=os.path.dirname(__file__) + "/data/emission-input/", type=str, help='path of emission inputs (= fcdoutput of sumo simulation)')
#     parser.add_argument('--output', default=os.path.dirname(__file__) + "/data/emission-output/", type=str, help='path of emission outputs')
    
#     global args 
#     args = parser.parse_args()

    # net needed to convert x, z to lat, long 
    # global net
    # net = sumolib.net.readNet(args.input + "kirchheim-street-names.net.xml")

    # global raw_pollution_data
    # raw_pollution_data = {}

# net = None
# args = None
# raw_pollution_data = None
# net = sumolib.net.readNet(os.path.dirname(__file__) + "/data/emission-input/kirchheim-street-names.net.xml")

# proj = net.getGeoProj()
# offset = net.getLocationOffset()
# if any(offset):
#     adjust = lambda x, y, _dx=offset[0], _dy=offset[1]: x + _dx, y + _dy
# else:
#     adjust = lambda *a: a

# def latlong(x, y, _proj=proj, _adjust=adjust):
#     return = _proj(*_adjust(x, y), inverse=True)

# def fast_iter(context, func):
#     for _, elem in context:
#         func(elem)
#         elem.clear()
#         while elem.getprevious() is not None:
#             del elem.getparent()[0]
#     del context

# def aggregate(
#     vehicle,
#     _fields=_fields, 
#     _empty=_empty,
#     _get=itemgetter(*_fields, 'x', 'y'),
#     _conv=latlong,
# ):
#   values, x, y = map(float, _get(vehicle.attrib))
#     lng, lat = _conv(x, y)
#     data = raw_pollution_data.setdefault(f"{.4f},{.4f}", _empty)
#     for f, v in zip(_fields, values):
#         data[f] += v

# async def parse_emissions():
#     # root = etree.parse(args.output + "test_emission_output.xml")
#     xml_file = os.path.dirname(__file__) + "/data/emission-output/emission_output.xml"
#     caqi_context = etree.iterparse(xml_file, tag="vehicle")
#     timer_start = timer()
#     fast_iter(caqi_context, aggregate)

# async def get_caqi_data():
#     # init_env()
#     global raw_pollution_data
#     raw_pollution_data = {}
#     raw_pollution_data = await parse_emissions()
#     air = caqi.calc_indices(raw_pollution_data)
#     # with open(args.output + "temp_caqi.json", "w") as json_file:
#     #     json.dump(air, json_file)
#     return air

# def aggregate(vehicle):
#     veh_id = int(vehicle.attrib["id"])
#     veh_noise = float(vehicle.attrib["noise"])
#     veh_type = vehicle.attrib["type"]
#     veh_eclass = vehicle.attrib["eclass"]
#     veh_co2 = float(vehicle.attrib["CO2"])
#     veh_co = float(vehicle.attrib["CO"])
#     veh_hc = float(vehicle.attrib["HC"])
#     veh_speed = float(vehicle.attrib["speed"])
#     veh_nox = float(vehicle.attrib["NOx"]) 
#     veh_pmx = float(vehicle.attrib["PMx"]) # mg/s
#     veh_fuel = float(vehicle.attrib["fuel"]) # ml/s
#     veh_x = float(vehicle.attrib["x"])
#     veh_y = float(vehicle.attrib["y"])
#     print(net.getLocationOffset())
#     # print("x: %s - y: %s" % (veh_x, veh_y))
#     # lng, lat = net.convertXY2LonLat(float(vehicle.attrib["x"]), float(vehicle.attrib["y"]))

#     # coordinate = str(round(lat, 4)) + "," + str(round(lng, 4))
#     coordinate = str(veh_x) + "," + str(veh_y)
#     if coordinate in raw_pollution_data:
#         raw_pollution_data[coordinate]["CO2"] += veh_co2
#         raw_pollution_data[coordinate]["NOX"] += veh_nox
#         raw_pollution_data[coordinate]["PMX"] += veh_pmx
#         raw_pollution_data[coordinate]["CO"] += veh_co
#     else:
#         raw_pollution_data[coordinate] = {}
#         raw_pollution_data[coordinate]["CO2"] = veh_co2
#         raw_pollution_data[coordinate]["NOX"] = veh_nox
#         raw_pollution_data[coordinate]["PMX"] = veh_pmx
#         raw_pollution_data[coordinate]["CO"] = veh_co

# async def parse_emissions():
#     # root = etree.parse(args.output + "test_emission_output.xml")
#     xml_file = args.output + "emission_output.xml"
#     caqi_context = etree.iterparse(xml_file, tag="vehicle")
#     # emission_cycle = {}
#     # emission_sum = {}
#     # max_co2 = 0
#     timer_start = timer()
#     fast_iter(caqi_context, aggregate)
#     seconds = timer() - timer_start

#     print("[etree] Finished parsing XML in %s seconds" % seconds)
#     return raw_pollution_data
    # print(raw_pollution_data)
    # root = etree.Element("root") 
    # for _, timestep in etree.iterparse(xml_file, tag="timestep"):
    #     step_id = float(timestep.attrib["time"])
    #     print("* step - %s -- %s" % (step_id, (step_id / 11800) * 100))
    #     # emission_cycle[step_id] = {}
    #     # emission_cycle[step_id]["step_id"] = step_id
    #     for vehicle in timestep.iterfind("vehicle"):
    #         veh_id = int(vehicle.attrib["id"])
    #         veh_noise = float(vehicle.attrib["noise"])
    #         veh_type = vehicle.attrib["type"]
    #         veh_eclass = vehicle.attrib["eclass"]
    #         veh_co2 = float(vehicle.attrib["CO2"])
    #         veh_co = float(vehicle.attrib["CO"])
    #         veh_hc = float(vehicle.attrib["HC"])
    #         veh_speed = float(vehicle.attrib["speed"])
    #         veh_nox = float(vehicle.attrib["NOx"]) 
    #         veh_pmx = float(vehicle.attrib["PMx"]) # mg/s
    #         veh_fuel = float(vehicle.attrib["fuel"]) # ml/s
    #         lng, lat = net.convertXY2LonLat(float(vehicle.attrib["x"]), float(vehicle.attrib["y"]))
            # print(veh_id)
            # print("** veh - %s" % veh_id)
            # print("*** - calculating coordinates")

            # coordinate = str(round(lat, 4)) + "," + str(round(lng, 4))
            # if coordinate in raw_pollution_data:
            #     raw_pollution_data[coordinate]["CO2"] += veh_co2
            #     raw_pollution_data[coordinate]["NOX"] += veh_nox
            #     raw_pollution_data[coordinate]["PMX"] += veh_pmx
            #     raw_pollution_data[coordinate]["CO"] = veh_co
            # else:
            #     raw_pollution_data[coordinate] = {}
            #     raw_pollution_data[coordinate]["CO2"] = veh_co2
            #     raw_pollution_data[coordinate]["NOX"] = veh_nox
            #     raw_pollution_data[coordinate]["PMX"] = veh_pmx
            #     raw_pollution_data[coordinate]["CO"] = veh_co
            
            # if raw_pollution_data[coordinate]["CO2"] > max_co2:
            #     max_co2 = raw_pollution_data[coordinate]["CO2"]
            
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
        # timestep.clear(keep_tail=True)
    # print(raw_pollution_data)
    # seconds = timer() - timer_start
    # print("[etree] Finished parsing XML in %s seconds" % seconds)
    # return raw_pollution_data
    # return caqi.calc_sub_indices(raw_pollution_data)
    # return calc_sub_indices(raw_pollution_data)
    # with open(args.output + "emission_cycle.json", "w") as json_file:
    #     json.dump(heatmap_data, json_file)
    # return "yes"
    # return heatmap_data

# async def get_caqi_data():
#     # parser = argparse.ArgumentParser()
#     # parser.add_argument('--output', default=os.path.dirname(__file__) + "/data/emission-output/", type=str, help='path of emission outputs')
#     # args = parser.parse_args()
#     init_env()
#     raw_pollution_data = await parse_emissions()
#     air = caqi.calc_indices(raw_pollution_data)
#     # with open(args.output + "temp_caqi.json", "w") as json_file:
#     #     json.dump(air, json_file)
#     return air

