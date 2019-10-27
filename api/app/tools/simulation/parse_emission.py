
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
# from app.db.mongodb import DB
from app.tools.simulation import calc_caqi as caqi
from app.tools.simulation import preprocessor as ip
from app.core.config import DEFAULT_NET_INPUT, EMISSION_OUTPUT_BASE
from app.crud.emissions import (
    get_caqi_emissions_for_sim,
    insert_caqi_emissions,
    insert_raw_emissions
)
from app.db.mongodb import AsyncIOMotorClient

if 'SUMO_HOME' in os.environ:
    tools = os.path.join(os.environ['SUMO_HOME'], 'tools')
    sys.path.append(tools)
else:   
    sys.exit("please declare environment variable 'SUMO_HOME'")

import sumolib
net = sumolib.net.readNet(DEFAULT_NET_INPUT)

class Parser():
    def __init__(self, db: AsyncIOMotorClient, simulation_id):
        self.db = db
        self.sim_id = simulation_id
        self.sim_output_path = EMISSION_OUTPUT_BASE + "emission_output_%s.xml" % self.sim_id

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

    async def get_caqi_data(self):
        caqi = await get_caqi_emissions_for_sim(self.db, self.sim_id)
        if caqi != None:
            print("[PARSER] Simulation has already been run. Fetching CAQI from DB...")
            return caqi["emissions"] 
        else: 
            print("[PARSER] parsing XML emission outputs from traffic simulation")
            timer_start = timer()
            emissions = self.parse_emissions(self.sim_output_path)
            seconds = timer() - timer_start
            print(emissions)
            print("[etree] Finished parsing XML in %s seconds" % seconds)
            
            print("[PARSER] Saving raw simulated emissions to database")
            raw_emissions = emissions.reset_index().to_json(orient='index')
            await insert_raw_emissions(self.db, self.sim_id, raw_emissions)

            print("[CAQI] calculating subindices and overall CAQI")
            caqi_emissions = emissions.apply(caqi.calc_indices, axis=1)
            result = caqi_emissions.reset_index().to_json(orient='index')

            print("[PARSER] Saving calculated inidzes to database")
            await insert_caqi_emissions(self.db, self.sim_id, result)
            return result