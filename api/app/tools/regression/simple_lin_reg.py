import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import datetime
import json

from fastapi import Depends
from app.db.mongodb import AsyncIOMotorClient, get_database
from app.crud.emissions import get_raw_emissions_from_sim
from app.core.config import (
    WEATHER_BASEDIR,
    WEATHER_PRESSURE,
    WEATHER_TEMP_HUMID, 
    WEATHER_WIND,
    AIR_BASEDIR,
    PLOT_BASEDIR
)

import numpy as np
from sklearn.linear_model import LinearRegression

class LinReg():
    def __init__(self, db: AsyncIOMotorClient, sim_id):
        self.db = db
        self.sim_id = sim_id
        self.raw_emission_columns = ['CO', 'NOx', 'PMx']
    
    async def fetch_simulated_emissions(self):
        # NOTE: First we fetch the simulated emissions
        raw_emissions = await get_raw_emissions_from_sim(self.db, self.sim_id)
        df = pd.DataFrame(pd.read_json(raw_emissions["emissions"], orient='index'))
        df_raw_em = df[self.raw_emission_columns]
        # NOTE: The simulation runs between 10000 - 20000 seconds so aggregate for every hour
        df_sim_nox = df_raw_em['NOx'].groupby(df.index // 3600).mean()
        df_sim_nox = pd.DataFrame({'NOx': df_sim_nox})
        # print(df_sim_nox)

        # NOTE: Fetch weather data and format timestamp
        df_temp_humid = self.get_temp_humid()
        df_temp = self.format_weather_by_key(
            df=df_temp_humid,
            key='TEMP',
            start_date='2019-02-01',
            end_date='2019-06-01',
            start_hour='06:00',
            end_hour='11:00'
        )
        # df_temp.plot(figsize=(18, 5))
        # plt.savefig(PLOT_BASEDIR + '/temp_02-06_morning')

        df_humidity = self.format_weather_by_key(
            df=df_temp_humid,
            key='HUMIDITY',
            start_date='2019-02-01',
            end_date='2019-06-01',
            start_hour='06:00',
            end_hour='11:00'
        )
        # df_humidity.plot(figsize=(18, 5))
        # plt.savefig(PLOT_BASEDIR + '/humidity_02-06_morning')

        df_pressure = self.get_pressure()
        df_pressure_nn = self.format_weather_by_key(
            df=df_pressure, 
            key='PRESSURE_NN', 
            start_date='2019-02-01', 
            end_date='2019-06-01', 
            start_hour='06:00', 
            end_hour='11:00'
        )
        # df_pressure_nn.plot(figsize=(18, 5))
        # plt.savefig(PLOT_BASEDIR + '/pressure-nn_02-06_morning')

        df_wind = self.get_wind()
        df_wind_speed = self.format_weather_by_key(
            df=df_wind, 
            key='WIND_SPEED', 
            start_date='2019-02-01', 
            end_date='2019-06-01', 
            start_hour='06:00', 
            end_hour='11:00'
        )
        
        # NOTE: Get real weather data and format it accordingly. Here we'll look at 2019 from 7:00 to 10:00
        air_df = self.get_real_air()
        # df_ex = self.format_real_air_by_key(
        #     df=air_df, 
        #     key='no2', 
        #     start_date='2019-02-15', 
        #     end_date='2019-02-16', 
        #     start_hour='0:00', 
        #     end_hour='23:00'
        # )
        # df_ex.plot(figsize=(18, 5))
        # plt.savefig(PLOT_BASEDIR + '/no2_15_02_whole_day')
        df_real_no2 = self.format_real_air_by_key(
            df=air_df, 
            key='no2', 
            start_date='2019-02-01', 
            end_date='2019-06-01', 
            start_hour='6:00', 
            end_hour='11:00'
        )
        # df_real_no2.plot(figsize=(18, 5))
        # plt.savefig(PLOT_BASEDIR + '/no2_02-06_morning')
        
        # print(df_real_no2)
        df_combined = pd.concat([df_temp, df_humidity, df_pressure_nn, df_real_no2], axis=1).dropna()
        length, _ = df_combined.shape
        df_sim_nox_repeated = pd.concat([df_sim_nox]*int(abs(length / 6)), ignore_index=True)
        df_combined = df_combined.assign(SIM_NOX=df_sim_nox_repeated.values)
        
        print(df_combined)
        # print(df_combined.shape)
        # df_combined.plot(figsize=(18, 5))
        # plt.savefig(PLOT_BASEDIR + '/combined-hum-temp-no2_sim-nox_02-06_morning')
        
        # NOTE: Linear regression
        X = df_combined[['TEMP', 'HUMIDITY', 'PRESSURE_NN', 'SIM_NOX']]
        Y = df_combined['no2']

        regr = LinearRegression()
        regr.fit(X, Y)

        print('Intercept: \n', regr.intercept_)
        print('Coefficients: \n', regr.coef_)







        df_temp_test = self.format_weather_by_key(
            df=df_temp_humid,
            key='TEMP',
            start_date='2019-09-01',
            end_date='2019-10-01',
            start_hour='06:00',
            end_hour='11:00'
        )
        # df_temp.plot(figsize=(18, 5))
        # plt.savefig(PLOT_BASEDIR + '/temp_02-06_morning')

        df_humidity_test = self.format_weather_by_key(
            df=df_temp_humid,
            key='HUMIDITY',
            start_date='2019-09-01',
            end_date='2019-10-01',
            start_hour='06:00',
            end_hour='11:00'
        )
        # df_humidity.plot(figsize=(18, 5))
        # plt.savefig(PLOT_BASEDIR + '/humidity_02-06_morning')

        df_pressure_test = self.format_weather_by_key(
            df=df_pressure, 
            key='PRESSURE_NN', 
            start_date='2019-09-01', 
            end_date='2019-10-01', 
            start_hour='06:00', 
            end_hour='11:00'
        )
        df_real_no2_test = self.format_real_air_by_key(
            df=air_df, 
            key='no2', 
            start_date='2019-09-01', 
            end_date='2019-10-01', 
            start_hour='6:00', 
            end_hour='11:00'
        )

        df_combined_test = pd.concat([df_temp_test, df_humidity_test, df_pressure_test, df_real_no2_test], axis=1).dropna()
        length_test, _ = df_combined_test.shape
        print('fucking %d' % length_test)
        print(abs(length_test / 6))
        df_sim_nox_repeated_test = pd.concat([df_sim_nox]*int(abs(length_test / 6)), ignore_index=True)
        df_combined_test = df_combined_test.assign(SIM_NOX=df_sim_nox_repeated_test.values)

        Z = df_combined_test[['TEMP', 'HUMIDITY', 'PRESSURE_NN', 'SIM_NOX']]
        print(regr.predict(Z))

        df_combined_test = df_combined_test.assign(predicted=regr.predict(Z))
        print(df_combined_test)
        df_combined_test[['predicted', 'no2']].plot(figsize=(18, 5))
        plt.savefig(PLOT_BASEDIR + '/prediction-09-10_morning')
        return raw_emissions

    ################################## WEATHER FUNCTIONS ########################################
    # weather inputs are: 
    # wind speed: column FF, in m/s
    # wind direction column D, degrees
    # pressure: column p0, to sea level (NN) reduced
    # relative humidty: column RF_TU, %
    # temperature: column TT_TU, °C
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

    def format_weather_by_key(self, df, key, start_date, end_date, start_hour, end_hour):
        df = df[['MESS_DATUM', key]]
        df['MESS_DATUM'] = df['MESS_DATUM'].apply(lambda x: datetime.datetime.strptime(str(x), '%Y%m%d%H'))
        # return df.set_index('MESS_DATUM')
        mask = (df['MESS_DATUM'] > start_date) & (df['MESS_DATUM'] <= end_date)
        df = df.loc[mask].set_index('MESS_DATUM')
        return df.between_time(start_hour, end_hour)
    
    #################################AIR FUNCTIONS#########################################
    def format_real_air_by_key(self, df, key, start_date, end_date, start_hour, end_hour):
        df = pd.DataFrame(df[key].tolist())
        df['time'] = pd.to_datetime(df['time'])
        df = df[['time', 'value']]
        mask = (df['time'] > start_date) & (df['time'] <= end_date)
        df = df.loc[mask].set_index('time')
        df = df.rename(columns={ 'value': key })
        return df.between_time(start_hour, end_hour)

    def get_real_air(self):
        air_frames = [self.get_real_air_from_file(AIR_BASEDIR + '/air_2019_%0*d.json' % (2, index)) for index in range(1, 11)]
        return pd.concat(air_frames, keys=['%d' % index for index in range(1, 11)]).dropna()
    
    def get_real_air_from_file(self, filepath):
        data = json.load(open(filepath))
        return pd.DataFrame(
            dict([ (k, pd.Series(v)) for k, v in data['features'][6]['properties']['timeValueSeries'].items() ])
        )

