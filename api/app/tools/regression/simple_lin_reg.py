import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import datetime
import json

from fastapi import Depends
from app.db.mongodb import AsyncIOMotorClient, get_database
from app.crud.emissions import (
    get_raw_emissions_from_sim, 
    fetch_air_traffic_from_hawa_dawa, 
    get_all_air_traffic
)
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
    def __init__(self, db: AsyncIOMotorClient, sim_id, existing_regr=None):
        self.db = db
        self.sim_id = sim_id
        self.existing_regr = existing_regr
        self.raw_emission_columns = ['CO', 'NOx', 'PMx']
        self.real_emission_columns = ['no2', 'pm2.5', 'pm10', 'o3']

    async def predict_emission(self, data=None, input_keys=['veh', 'TEMP', 'HUMIDITY', 'PMx'], output_key='pm10'):
        df_combined = await self.aggregate_data('2019-08-22', '2019-10-15', '6:00', '11:00')
        df_combined = df_combined.dropna()
        # print(df_combined)
        regr = await self.train_model(df_combined, input_keys, output_key)

        df_test = await self.aggregate_data('2019-10-16', '2019-10-22', '6:00', '11:00')
        df_test = df_test.dropna()
        Z = df_test[input_keys]
        # print(regr.predict(Z))
        # df_test = df_test.assign({output_key + '_predicted': regr.predict(Z)})
        df_test[output_key + '_predicted'] = regr.predict(Z)
        print(df_test)
        # self.save_df_to_plot(df_test[[output_key, 'predicted']], '%s_prediction_new' % output_key)
        df_test = df_test.reset_index()
        return df_test[[output_key, '%s_predicted' % output_key]]

    async def train_model(self, df, input_keys, output_key):
        X = df[input_keys]
        Y = df[output_key]

        regr = LinearRegression()
        regr.fit(X, Y)

        print('Intercept: \n', regr.intercept_)
        print('Coefficients: \n', regr.coef_)
        return regr



    async def get_hw_data(self):
        # data = await fetch_air_traffic_from_hawa_dawa(self.db, start_date='2019-10-01', end_date='2019-10-20')
        # result = await get_all_air_traffic(self.db)
        # df = pd.DataFrame(result)
        data = await self.get_real_air()
        print(data)
        # print(df)

    ######################################################################################################
    ################################## DATA AGGREGATION FUNCTIONS ########################################
    ######################################################################################################
    async def aggregate_data(self, start_date='2019-08-16', end_date='2019-10-24', start_hour='6:00', end_hour='11:00'):
        df_sim = await self.fetch_simulated_emissions()
        df_weather = await self.fetch_weather_data(start_date, end_date, start_hour, end_hour)
        df_air_traffic = await self.fetch_air_and_traffic(start_date, end_date, start_hour, end_hour)
        # print(df_air_traffic)
        df_combined = pd.concat([df_air_traffic, df_weather], axis=1)
        # print(df_combined)

        # NOTE: Combine all dataframes
        length, _ = df_combined.shape
        df_sim_repeated = pd.concat([df_sim]*int(abs((length / 6))), ignore_index=True)
        df_combined = df_combined.reset_index().rename(columns={'index': 'time'}) # prepare for concat with sim values
        return pd.concat([df_combined, df_sim_repeated], axis=1).set_index(['time'])


    #####################################################################################################
    ################################## DATA COLLECTION FUNCTIONS ########################################
    ######################################################################################################
    async def fetch_simulated_emissions(self):
        # NOTE: First we fetch the simulated emissions
        raw_emissions = await get_raw_emissions_from_sim(self.db, self.sim_id)
        df = pd.DataFrame(pd.read_json(raw_emissions["emissions"], orient='index'))
        df_raw_em = df[self.raw_emission_columns]
        # NOTE: The simulation runs between 10000 - 20000 seconds so aggregate for every hour
        # df_sim_nox = df_raw_em['NOx'].groupby(df.index // 3600).mean()
        return df_raw_em.groupby(df.index // 3600).mean()

    async def fetch_weather_data(self, start_date='2019-01-01', end_date='2019-10-28', start_hour='0:00', end_hour='23:00'):
        # NOTE: Fetch weather data and format timestamp
        df_temp_humid = self.get_temp_humid()
        df_temp = self.format_weather_by_key(
            df=df_temp_humid,
            key='TEMP',
            start_date=start_date,
            end_date=end_date,
            start_hour=start_hour,
            end_hour=end_hour
        )
        # self.save_df_to_plot(df_temp, 'temp_02-06_morning')

        df_humidity = self.format_weather_by_key(
            df=df_temp_humid,
            key='HUMIDITY',
            start_date=start_date,
            end_date=end_date,
            start_hour=start_hour,
            end_hour=end_hour
        )
        # self.save_df_to_plot(df_humidity, 'humidity_02-06_morning')

        df_pressure = self.get_pressure()
        df_pressure_nn = self.format_weather_by_key(
            df=df_pressure, 
            key='PRESSURE_NN', 
            start_date=start_date, 
            end_date=end_date, 
            start_hour=start_hour, 
            end_hour=end_hour
        )
        # self.save_df_to_plot(df_pressure_nn, 'pressure-nn_02-06_morning')

        df_wind = self.get_wind()
        df_wind_speed = self.format_weather_by_key(
            df=df_wind, 
            key='WIND_SPEED', 
            start_date=start_date, 
            end_date=end_date, 
            start_hour=start_hour, 
            end_hour=end_hour
        )
        # self.save_df_to_plot(df_wind_speed, 'pressure-nn_02-06_morning')
        return pd.concat([frame for frame in [df_temp, df_humidity, df_pressure_nn, df_wind_speed] if not frame.empty], axis=1)
        

    async def fetch_air_and_traffic(self, start_date='2019-01-01', end_date='2019-10-28', start_hour='0:00', end_hour='23:00'):
        # NOTE: Get real weather data and format it accordingly. Here we'll look at 2019 from 7:00 to 10:00
        air_df = await self.get_real_air()

        bremicker_df = await self.get_bremicker()
        df_traffic = await self.format_real_air_by_key(
            bremicker_df,
            'veh',
            start_date, 
            end_date, 
            start_hour, 
            end_hour
        )
        # self.save_df_to_plot(df_traffic, 'bremicker_all')
        # print(df_traffic)
        frames = [
            await self.format_real_air_by_key(
                air_df,
                key,
                start_date, 
                end_date, 
                start_hour, 
                end_hour) for key in self.real_emission_columns
        ]

        if not df_traffic.empty:
            frames.append(df_traffic)

        return pd.concat(frames, axis=1)

    def save_df_to_plot(self, df, filename):
        if not df.empty:
            df.plot(figsize=(18, 5))
            plt.savefig(PLOT_BASEDIR + '/' + filename)
        else:
            print('[PLOT] Error saving plot. Dataframe empty!')

    #############################################################################################
    ################################## WEATHER FUNCTIONS ########################################
    #############################################################################################
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
    
    #######################################################################################
    #################################AIR FUNCTIONS#########################################
    #######################################################################################
    async def format_real_air_by_key(self, df, key, start_date, end_date, start_hour, end_hour):
        df = pd.DataFrame(df[key].tolist())
        df['time'] = pd.to_datetime(df['time'])
        df = df[['time', 'value']]
        mask = (df['time'] > start_date) & (df['time'] <= end_date)
        df = df.loc[mask].set_index('time')
        df = df.rename(columns={ 'value': key })
        return df.between_time(start_hour, end_hour)

    async def get_real_air(self):
        air_frames = [await self.get_real_air_from_file(AIR_BASEDIR + '/air_2019_%0*d.json' % (2, index)) for index in range(1, 11)]
        return pd.concat(air_frames, keys=['%d' % index for index in range(1, 11)]).dropna()
        # air_frames = await get_all_air_traffic(self.db)
        # data = [pd.DataFrame.from_csv(frame['data']) for frame in air_frames]
        # return data
        # for frame in air_frames:
        #     data = json.load(frame['data'])
        #     print(data)
        # 
    
    async def get_real_air_from_file(self, filepath):
        data = json.load(open(filepath))
        return pd.DataFrame(
            dict([ (k, pd.Series(v)) for k, v in data['features'][6]['properties']['timeValueSeries'].items() ])
        )
    
    ###############################################################################################
    ################################## BREMICKER FUNCTIONS ########################################
    ###############################################################################################
    async def get_bremicker(self):
        traffic_frames = [await self.get_bremicker_from_file(AIR_BASEDIR + '/air_2019_%0*d.json' % (2, index)) for index in range(1, 11)]
        return pd.concat(traffic_frames, keys=['%d' % index for index in range(1, 11)]).dropna()
    
    async def get_bremicker_from_file(self, filepath):
        data = json.load(open(filepath))
        return pd.DataFrame(
            dict([ (k, pd.Series(v)) for k, v in data['features'][2]['properties']['timeValueSeries'].items() ])
        )