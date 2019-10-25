
import os
import sys
import json
import asyncio
import pandas as pd
import numpy as np
import datetime

from lxml import etree
from operator import itemgetter
from timeit import default_timer as timer
from database.database import DB
from simulation import calc_caqi as caqi
from simulation import preprocessor as ip
from simulation.constants import Constants

if 'SUMO_HOME' in os.environ:
    tools = os.path.join(os.environ['SUMO_HOME'], 'tools')
    sys.path.append(tools)
else:   
    sys.exit("please declare environment variable 'SUMO_HOME'")

import sumolib
net = sumolib.net.readNet(Constants.DEFAULT_NET_INPUT)

class Parser():
    def __init__(self, simulation_id):
        self.sim_id = simulation_id
        self.sim_output_path = Constants.EMISSION_OUTPUT_BASE + "emission_output_%s.xml" % self.sim_id

    def extract_attributes(self, context, fields):
        values = itemgetter(*fields)
        for _, elem in context:
            yield values(elem.attrib)
            elem.clear()
            while elem.getprevious() is not None:
                del elem.getparent()[0]
        del context

    def parse_emissions(self, filepath):
        context = etree.iterparse(filepath, tag="vehicle")

        # create a dataframe from XML data a single call
        coords = ['x', 'y']
        entries = ['CO2', 'CO', 'NOx', 'PMx', 'fuel']
        df = pd.DataFrame(
            self.extract_attributes(context, coords + entries),
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

    def get_caqi_data(self):
        print("[PARSER] parsing XML emission outputs from traffic simulation")
        timer_start = timer()
        emissions = self.parse_emissions(self.sim_output_path)
        seconds = timer() - timer_start
        print(emissions)
        print("[etree] Finished parsing XML in %s seconds" % seconds)
        
        print("[CAQI] calculating subindices and overall CAQI")
        caqi_emissions = emissions.apply(caqi.calc_indices, axis=1)
        result = caqi_emissions.reset_index().to_json(orient='index')
        print("[PARSER] Saving calculated inidzes to database")
        DB.insert(collection='caqi_emissions', data={ "created_at:": datetime.datetime.utcnow(), "sim_id": self.sim_id, "emissions": result })
        return result
        # return result.to_json(orient='index')
        # return json.dumps(result)
        # print(json_val)
        # return caqi_emissions.to_json(orient='index')
        # df = pd.DataFrame(emissions)
        # df.app
        # return emissions