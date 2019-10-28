import numpy as np 
import pandas as pd
import json

from fastapi import Depends
from app.db.mongodb import AsyncIOMotorClient, get_database
from app.crud.emissions import get_raw_emissions_from_sim
from app.core.config import (
    WEATHER_BASEDIR,
    WEATHER_PRESSURE,
    WEATHER_TEMP_HUMID, 
    WEATHER_WIND,
    AIR_BASEDIR
)

import numpy as np
from sklearn.linear_model import LinearRegression

class LinReg():
    def __init__(self, db: AsyncIOMotorClient, sim_id):
        self.db = db
        self.sim_id = sim_id
        self.raw_emission_columns = ['CO', 'NOx', 'PMx']
    
    async def fetch_simulated_emissions(self):
        raw_emissions = await get_raw_emissions_from_sim(self.db, self.sim_id)
        # print(raw_emissions["emissions"])
        pd_emissions = pd.read_json(raw_emissions["emissions"], orient='index')
        # df = pd.DataFrame(pd_emissions, dtype=float)
        df = pd.DataFrame(pd_emissions)
        # print(df[self.raw_emission_columns])

        df_temp_humid = self.get_temp_humid()
        df_pressure = self.get_pressure()
        df_wind = self.get_wind()
        df_air = self.get_real_air(AIR_BASEDIR + '/air_2019_01.json')
        return raw_emissions

    # weather inputs are: 
    # wind speed: column FF, in m/s
    # wind direction column D, degrees
    # pressure: column p0, to sea level (NN) reduced
    # relative humidty: column RF_TU, %
    # temperature: column TT_TU, °C

    # def get_weather_as_np(self):
        # with open()

    def get_temp_humid(self):
        df = pd.read_csv(WEATHER_TEMP_HUMID, delimiter=';')
        df.columns = df.columns.str.strip()
        df = df.rename(columns={"TT_TU": "TEMP", "RF_TU": "HUMIDITY"})
        return df[['MESS_DATUM', 'TEMP', 'HUMIDITY']]
        # humidty, temperature = np.loadtxt(WEATHER_TEMP_HUMID, delimiter=';', usecols=(3, 4), skiprows=1, unpack=True)

    def get_pressure(self):
        df = pd.read_csv(WEATHER_PRESSURE, delimiter=';')
        df.columns = df.columns.str.strip()
        df = df.rename(columns={"P": "PRESSURE_NN", "P0": "PRESSURE_STATION"})
        return df[['MESS_DATUM', 'PRESSURE_NN', 'PRESSURE_STATION']]
        # return 0

    def get_wind(self):
        df = pd.read_csv(WEATHER_WIND, delimiter=';')
        df.columns = df.columns.str.strip()
        df = df.rename(columns={"FF": "WIND_SPEED", "DD": "WIND_DIR"})
        return df[['MESS_DATUM', 'WIND_SPEED', 'WIND_DIR']]
    

    def get_real_air(self, filepath):
        data = json.load(open(filepath))

        df = pd.DataFrame(
            dict([ (k, pd.Series(v)) for k,v in data['features'][6]['properties']['timeValueSeries'].items() ])
        )
        return df

