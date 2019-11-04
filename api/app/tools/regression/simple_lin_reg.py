import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import datetime
import json

from fastapi import Depends
from app.db.mongodb import AsyncIOMotorClient, get_database
from app.crud.emissions import (
    get_raw_emissions_from_sim
)
from app.crud.hawa_dawa import (
    # fetch_data_by_month_from_hawa_dawa, 
    get_hawa_dawa_by_time,
    # get_all_hawa_dawa
)
from app.tools.regression.utils.weather import (
    fetch_weather_data
)
from app.tools.simulation.parse_emission import Parser
from app.tools.regression.utils.bremicker_boxes import bremicker_boxes
import app.tools.simulation.calc_caqi as aqi
from app.core.config import (
    WEATHER_BASEDIR,
    WEATHER_PRESSURE,
    WEATHER_TEMP_HUMID, 
    WEATHER_WIND,
    AIR_BASEDIR,
    PLOT_BASEDIR,
    EMISSION_OUTPUT
)

import numpy as np
from sklearn.linear_model import LinearRegression

class LinReg():
    def __init__(self, db: AsyncIOMotorClient, sim_id=None, existing_regr=None):
        self.db = db
        self.sim_id = sim_id
        self.existing_regr = existing_regr
        self.raw_emission_columns = ['CO', 'NOx', 'PMx']
        self.real_emission_columns = ['no2', 'pm2.5', 'pm10', 'o3']


    # async def get_hw_data(self):
    #     traffic = await get_all_air_traffic(self.db)
    #     print(traffic)
        # for element in traffic:
        #     df_t = pd.DataFrame(element)
        #     print(df_t)
        # df = pd.DataFrame(traffic)
        # df_traffic = pd.DataFrame(df['data'])
        # print(df_traffic)
    async def get_sim_em_distribution(self):
        
        df_sim = await self.fetch_simulated_emissions(672)
        # print(df_sim.shape)
        # print(df_sim)
        # return df_sim
        # .apply(aqi.calc_nox_index)
        # df_nox = df_sim['NOx']
        # print(df_nox)
        # var = df_nox.var(axis=0)
        # std = df_nox.std(axis=0)
        # print(var)
        # print(std)
        # print(df_nox)
    # , figsize=(18, 5)
        # df_nox.plot.hist(grid=True, bins=300, rwidth=0.9, color='#333333')
        # plt.savefig(PLOT_BASEDIR + '/' + 'hist_nox')


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

    ######################################################################################################
    ################################## DATA AGGREGATION FUNCTIONS ########################################
    ######################################################################################################
    async def aggregate_data(self, start_date='2019-08-16', end_date='2019-10-24', start_hour='6:00', end_hour='11:00'):
        df_sim = await self.fetch_simulated_emissions(672)
        df_weather = await fetch_weather_data(start_date, end_date, start_hour, end_hour)
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
    async def fetch_simulated_emissions(self, box_id):
        # NOTE: First we fetch the simulated emissions
        parser = Parser(self.db, self.sim_id)
        lat, lng = [round(bremicker_boxes[box_id]['lat'], 3), round(bremicker_boxes[box_id]['lng'], 3)]
        # raw_emissions = await get_raw_emissions_from_sim(self.db, self.sim_id)
        raw_emissions = await parser.get_caqi_data()
        print(raw_emissions)
        # df = pd.DataFrame(pd.read_json(raw_emissions["emissions"], orient='index'))
        # print(df)
        # print(lat, lng)
        # latlng = ['lat', 'lng']
        # df = df[latlng + self.raw_emission_columns]
        # print(df)
        # df = df.groupby([df.index // 60] + latlng)[self.raw_emission_columns].mean().reset_index(latlng)
        
        # mask = (round(df['lat'], 3) == lat) & (round(df['lng'], 3) == lng)
        # df = df.loc[mask]
        # return df
        # print(df)
        # df[latlng] = df[latlng].round(3)
        # df = df.groupby([df.index] + latlng)[self.raw_emission_columns].sum()
        # print(df)
        # df[latlng] = df[latlng].round(3)

        # aggregate the results and return summed dataframe
        # df = df.groupby(latlng)[self.raw_emission_columns].sum()
        # print(df)
        # return df_raw_em.to_dict(orient='index')


        # df_raw_em = df_raw_em[df_raw_em.lat == lat, df_raw_em.lng == lng]
        # print(df_raw_em)
        # NOTE: The simulation runs between 10000 - 20000 seconds so aggregate for every hour
        # df_sim_nox = df_raw_em['NOx'].groupby(df.index // 3600).mean()
        # return df_raw_em.groupby(df.index // 60)[self.raw_emission_columns].mean()


        

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

    async def get_mean_vehicle_by_hour(self, start_date, end_date, start_hour, end_hour):
        bremicker_df = await self.get_bremicker()
        df_traffic = await self.format_real_air_by_key(
            bremicker_df,
            'veh',
            start_date, 
            end_date, 
            start_hour, 
            end_hour
        )
        # self.save_df_to_plot(df_traffic, 'big_number_vehicles_7-11')
        df_sum = df_traffic.reset_index()
        # print(df_avg)
        df_sum = df_sum.groupby(by=df_sum['time'].dt.date).sum()
        print(df_sum)
        print(df_sum.mean())
        return df_sum
        # self.save_df_to_plot(df_traffic, 'number_vehicles_sum_diesel_7am-11am')
        # print(df_sum)

        # print(df_traffic.mean())

    #############################################################################################
    ################################## WEATHER FUNCTIONS ########################################
    #############################################################################################
    # weather inputs are: 
    # wind speed: column FF, in m/s
    # wind direction column D, degrees
    # pressure: column p0, to sea level (NN) reduced
    # relative humidty: column RF_TU, %
    # temperature: column TT_TU, °C

    
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
            dict([ (k, pd.Series(v)) for k, v in data['features'][1]['properties']['timeValueSeries'].items() ])
        )
    
    async def get_bremicker_sensors_from_file(self, filepath):
        data = json.load(open(filepath))
        traffic_frames = []
        for feature in data['features']:
            if feature['properties']['type'] == 'bremicker':
                sensor = {}
                sensor['coordinates'] = feature['geometry']['coordinates']
                df = pd.DataFrame(
                    dict([ (k, pd.Series(v)) for k, v in feature['properties']['timeValueSeries'].items() ])
                )
                sensor['vehicleNumber'] = feature['properties']['timeValueSeries'].items()
                traffic_frames.append(sensor)
        return traffic_frames




    async def get_air_sensor_data(self):
        # result = await get_all_hawa_dawa(self.db)
        result = await get_hawa_dawa_by_time(self.db)
        print(result)
        # months = []
        # for elem in result:
        #     data = json.loads(elem['data'])
        #     for feature in data['features']:
        #         if feature['properties']['type'] == 'hawadawa' and feature['properties']['timeValueSeries']:
        #             df = pd.DataFrame(dict([ (k, pd.Series(v)) for k, v in feature['properties']['timeValueSeries'].items() ]))
        #             # df = df.apply(self.format_column, result_type='broadcast')
        #             frames = []
        #             for pollutant in df.columns:
        #                 df_pol = df[pollutant].apply(pd.Series)[['time', 'value']].dropna()
        #                 df_pol['time'] = pd.to_datetime(df_pol['time'])
        #                 df_pol = df_pol.set_index('time')
        #                 df_pol = df_pol.rename(columns={'value': pollutant})
        #                 frames.append(df_pol)
        #             months.append( pd.concat(frames, axis=1) )
        # result = pd.concat(months)
        # print(result)
        # return result
    # def format_column(self, col):
    #     df_pol = col.apply(pd.Series)[['time', 'value']].dropna()
    #     df_pol['time'] = pd.to_datetime(df_pol['time'])
    #     df_pol = df_pol.set_index('time')
    #     df_pol = df_pol.rename(columns={'value': col.name})
    #     # print(df_pol)
    #     return df_pol
        # result = []
        # for df in frames:
        #     frames = []
        #     for pollutant in df.columns:
        #         df_pol = df[pollutant].apply(pd.Series)[['time', 'value']].dropna()
        #         df_pol['time'] = pd.to_datetime(df_pol['time'])
        #         df_pol = df_pol.set_index('time')
        #         df_pol = df_pol.rename(columns={'value': pollutant})
        #         frames.append(df_pol)
        #     result.append( pd.concat(frames, axis=1) )
        # result = pd.concat(result)
        # print(result)

        # result = await get_all_air_traffic(self.db)
        # frames = []
        # for elem in result:
        #     data = json.loads(elem['data'])
        #     for feature in data['features']:
        #         if feature['properties']['type'] == 'hawadawa' and feature['properties']['timeValueSeries']:
        #             frames.append(
        #                 pd.DataFrame(dict([ (k, pd.Series(v)) for k, v in feature['properties']['timeValueSeries'].items() ]))
        #                 # pd.DataFrame(feature['properties']['timeValueSeries']).set_index('time')
        #             )
        # result = []
        # for df in frames:
        #     frames = []
        #     for pollutant in df.columns:
        #         df_pol = df[pollutant].apply(pd.Series)[['time', 'value']].dropna()
        #         df_pol['time'] = pd.to_datetime(df_pol['time'])
        #         df_pol = df_pol.set_index('time')
        #         df_pol = df_pol.rename(columns={'value': pollutant})
        #         frames.append(df_pol)
        #     result.append( pd.concat(frames, axis=1) )
        # result = pd.concat(result)
        # print(result)