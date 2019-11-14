import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import datetime
import json
import math
from tbats import BATS, TBATS

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
from app.crud.bremicker import (
    get_bremicker_by_time
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

from app.tools.regression.utils.weather import fetch_weather_data

import numpy as np
# from tensorflow import keras
from keras.models import Sequential
from keras.layers import Dense
from keras.layers import LSTM
from sklearn.linear_model import LinearRegression
from sklearn.metrics import (mean_absolute_error, mean_squared_error)
from sklearn.neural_network import MLPRegressor
from sklearn import preprocessing as pre

class LinReg():
    def __init__(self, db: AsyncIOMotorClient, sim_id=None, existing_regr=None):
        self.db = db
        self.sim_id = sim_id
        self.existing_regr = existing_regr
        self.raw_emission_columns = ['CO', 'NOx', 'PMx']
        self.real_emission_columns = ['no2', 'pm2.5', 'pm10', 'o3']

    
    async def start_cnn(
        self,
        start_date='2019-08-01', 
        end_date='2019-10-20', 
        start_hour='7:00', 
        end_hour='10:00', 
        data=None, 
        boxID=672, 
        input_keys=['temp', 'hum', 'PMx', 'WIND_SPEED', 'WIND_DIR'], 
        output_key='pm10'
    ):
        input_keys.append(boxID)
        df = await self.aggregate_data(boxID, start_date, end_date, start_hour, end_hour)
        # print('aggregated')
        # print(df)

        rows = round(df.shape[0] * 0.8)
        df_train = df.iloc[:rows]
        df_test = df.iloc[rows:]
        scaler = pre.StandardScaler()
        train_scaled = scaler.fit_transform(df.iloc[:rows][input_keys])
        test_scaled = scaler.fit_transform(df.iloc[rows:][input_keys])

        model = MLPRegressor(
            hidden_layer_sizes=(10,),
            activation='relu',
            solver='adam',
            learning_rate='adaptive',
            max_iter=1000,
            learning_rate_init=0.01,
            alpha=0.01
        )
        model.fit(train_scaled, df_train[output_key])
        
        # df_test[output_key + '_mlp_predicted'] = model.predict(df_test[input_keys])
        df_test[output_key + '_mlp_predicted'] = model.predict(test_scaled)
        result = df_test[[output_key, '%s_mlp_predicted' % output_key]]
        # print(result)
        print("Mean Abs Error MLP: " + str(mean_absolute_error(result[output_key].to_numpy(), result['%s_mlp_predicted' % output_key].to_numpy())))
        # self.save_df_to_plot(result, '%s_mlp_dist_regressor' % output_key.replace('.', '-'))
        return result


    def normalize_data(self, column):
        std = column.std(axis=0)
        var = column.var(axis=0)
        mean = column.mean(axis=0)
        print('Var: \n', var)
        print('Std: \n', std)
        print('Mean: \n', mean)
        return (column - mean) / std
    
    async def start_tbats(
        self,
        start_date='2019-08-01', 
        end_date='2019-10-20', 
        start_hour='7:00', 
        end_hour='10:00', 
        data=None, 
        boxID=672, 
        input_keys=['temp', 'hum', 'PMx', 'WIND_SPEED', 'WIND_DIR'], 
        output_key='pm10'
    ):
        input_keys.append(boxID)
        df = await self.aggregate_data(boxID, start_date, end_date, start_hour, end_hour)
        # print('aggregated')
        # print(df)

        rows = round(df.shape[0] * 0.8)
        df_train = df.iloc[:rows]
        df_test = df.iloc[rows:]
        # scaler = pre.StandardScaler()
        scaler = pre.MinMaxScaler(feature_range=(0, 100))
        df_train[['pm10']] = scaler.fit_transform(df_train[['pm10']])
        print(df_train)
        # test_scaled = scaler.fit_transform(df_test[input_keys])
        #[val for val in np.arange(4, 320, 80)]
        
        # periods = [val for val in np.arange(4, 256, 4)]
        # print(periods)
        # estimator = TBATS(seasonal_periods=periods)
        
        # # fitted_model = estimator.fit(df_train[['temp', 'hum', 'WIND_SPEED', 'WIND_DIR', 'pm10']])
        # fitted_model = estimator.fit(df_train[['pm10']])
        # y_forecasted = fitted_model.forecast(steps=64)
        # print(fitted_model.summary())

        # # Time series analysis
        # print(fitted_model.y_hat) # in sample prediction
        # print(fitted_model.resid) # in sample residuals
        # print(fitted_model.aic)

        # print(fitted_model.params.alpha)
        # print(fitted_model.params.beta)
        # print(fitted_model.params.x0)
        # print(fitted_model.params.components.use_box_cox)
        # print(fitted_model.params.components.seasonal_harmonics)
        # print(y_forecasted)
        # # print(df_test)
        # df_test['predicted'] = y_forecasted
        # print(df_test)
        # self.save_df_to_plot(df_test[['pm10', 'predicted']], 'tbats_pm10_forecast')
    # convert an array of values into a dataset matrix

    def create_dataset(self, dataset, look_back=1):
        dataX, dataY = [], []
        for i in range(len(dataset)-look_back-1):
            a = dataset[i:(i+look_back), 0]
            dataX.append(a)
            dataY.append(dataset[i + look_back, 0])
        return np.array(dataX), np.array(dataY)

    async def start_lstm(
        self,
        start_date='2019-08-01', 
        end_date='2019-10-20', 
        start_hour='7:00', 
        end_hour='10:00', 
        data=None, 
        boxID=672, 
        input_keys=['temp', 'hum', 'PMx', 'WIND_SPEED', 'WIND_DIR'], 
        output_key='pm10'
    ):
        input_keys.append(boxID)
        df = await self.aggregate_data(boxID, start_date, end_date, start_hour, end_hour)
        # feature_range=(0, 100)
        scaler = pre.MinMaxScaler()
        df_scaled = df
        df_scaled[[output_key]] = scaler.fit_transform(df[[output_key]])
        dataset = df_scaled[[output_key]]
        rows = round(df.shape[0] * 0.8)
        df_train = df_scaled[[output_key]].iloc[:rows]
        df_test = df_scaled[[output_key]].iloc[rows:]

        print(dataset)
        # self.save_df_to_plot(dataset, 'unscaled_%s' % output_key)
        # df_train = df_train.reset_index()
        # df_test = df_test.reset_index()
        # print(df_train)
        # print(df_train[['pm10']])

        look_back = 1
        trainX, trainY = self.create_dataset(df_train[[output_key]].values, look_back)
        testX, testY = self.create_dataset(df_test[[output_key]].values, look_back)

        # reshape input to be [samples, time steps, features]
        trainX = np.reshape(trainX, (trainX.shape[0], 1, trainX.shape[1]))
        testX = np.reshape(testX, (testX.shape[0], 1, testX.shape[1]))
        
        # create and fit the LSTM network
        model = Sequential()
        model.add(LSTM(4, input_shape=(1, look_back)))
        model.add(Dense(1))
        model.compile(loss='mean_squared_error', optimizer='adam')
        model.fit(trainX, trainY, epochs=100, batch_size=1, verbose=2)

        # make predictions
        trainPredict = model.predict(trainX)
        testPredict = model.predict(testX)
        # invert predictions
        trainPredict = scaler.inverse_transform(trainPredict)
        trainY = scaler.inverse_transform([trainY])
        testPredict = scaler.inverse_transform(testPredict)
        testY = scaler.inverse_transform([testY])
        # calculate root mean squared error
        trainScore = math.sqrt(mean_squared_error(trainY[0], trainPredict[:,0]))
        print('Train Score: %.2f RMSE' % (trainScore))
        testScore = math.sqrt(mean_squared_error(testY[0], testPredict[:,0]))
        print('Test Score: %.2f RMSE' % (testScore))

        # shift train predictions for plotting
        trainPredictPlot = np.empty_like(dataset)
        trainPredictPlot[:, :] = np.nan
        trainPredictPlot[look_back:len(trainPredict)+look_back, :] = trainPredict
        # shift test predictions for plotting
        testPredictPlot = np.empty_like(dataset)
        testPredictPlot[:, :] = np.nan
        testPredictPlot[len(trainPredict)+(look_back*2)+1:len(dataset)-1, :] = testPredict
        # plot baseline and predictions
        # print(scaler.inverse_transform(dataset))
        plt_dataset, = plt.plot(scaler.inverse_transform(dataset))
        plt_train_pred, = plt.plot(trainPredictPlot)
        plt_test_pred, = plt.plot(testPredictPlot)
        plt.legend([plt_dataset, plt_train_pred, plt_test_pred], ['%s Real Data' % output_key, '%s Training Prediction' % output_key, '%s Testing Prediction' % output_key])
        plt.show()


    # async def predict_emission(self, data=None, input_keys=['veh', 'TEMP', 'HUMIDITY', 'PMx'], output_key='pm10'):
    async def predict_emission(
        self,
        start_date='2019-08-01', 
        end_date='2019-10-20', 
        start_hour='7:00', 
        end_hour='10:00',
        data=None, 
        boxID=672, 
        input_keys=['temp', 'hum', 'PMx', 'WIND_SPEED', 'WIND_DIR'], 
        output_key='pm10'
    ):
        input_keys.append(boxID)

        df_combined = await self.aggregate_data(boxID, start_date, end_date, start_hour, end_hour)
        rows = round(df_combined.shape[0] * 0.8)
        df_train = df_combined.iloc[:rows]
        df_test = df_combined.iloc[rows:]
        
        # regr = await self.train_model(df_train, input_keys, output_key)

        if data != None:
            df_test = data
        
        model = LinearRegression()
        model.fit(df_train[input_keys], df_train[output_key])

        print('Intercept: \n', model.intercept_)
        print('Coefficients: \n', model.coef_)
        df_test[output_key + '_lin_predicted'] = model.predict(df_test[input_keys])
        # print(df_test)
        # self.save_df_to_plot(df_test[[output_key, '%s_predicted' % output_key]], '%s_lin_dist_prediction' % output_key)
        # df_test = df_test.reset_index()
        result = df_test[[output_key, '%s_lin_predicted' % output_key]]
        print("Mean Abs Error LinReg: " + str(mean_absolute_error(result[output_key].to_numpy(), result['%s_lin_predicted' % output_key].to_numpy())))
        print(result)
        return result

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
    async def aggregate_data(self, boxID=672, start_date='2019-08-01', end_date='2019-10-20', start_hour='7:00', end_hour='10:00'):
        df_air = await self.fetch_air_and_traffic(
            boxID, start_date, end_date, start_hour, end_hour
        )

        df_air.index.name = 'time'
        
        df_sim = await self.fetch_simulated_emissions(boxID)
        df_sim = df_sim.groupby('time')[self.raw_emission_columns].sum()
        # print(df_sim)

        df_nox = df_sim['NOx']
        df_pmx = df_sim['PMx']
        nox_var = df_nox.var(axis=0)
        nox_std = df_nox.std(axis=0)
        nox_mean = df_nox.mean(axis=0)
        pmx_var = df_pmx.var(axis=0)
        pmx_std = df_pmx.std(axis=0)
        pmx_mean = df_pmx.mean(axis=0)
        print('Var: \n', nox_var, pmx_var)
        print('Std: \n', nox_std, pmx_std)
        print('Mean: \n', nox_mean, pmx_mean)
        # print(df_air.shape)
        # ratio = ((nox_std / nox_mean) + (pmx_std / pmx_mean)) / 2
        sim_rows = df_sim.shape[0]
        air_rows = df_air.shape[0]
        rows_needed = (air_rows / (sim_rows / 60))
        ratio_nox = 1 - (nox_std / nox_mean)
        ratio_pmx = 1 - (pmx_std / pmx_mean)
        
        # print("%s rows needed" % rows_needed)
        # print("%s ratio_nox" % ratio_nox)
        # print("%s ratio_pmx" % ratio_pmx)
        
        nox_frames = [df_sim[['NOx']] * val for val in np.arange(1 - ratio_nox, 1 + ratio_nox, (ratio_nox * 2) / rows_needed)]
        df_nox = pd.concat(nox_frames, axis=0, ignore_index=True)
        df_nox = df_nox.groupby(df_nox.index // 60).sum()
        # self.save_df_to_plot(df_nox, 'nox_simulated_oversampled')

        pmx_frames = [df_sim[['PMx']] * val for val in np.arange(1 - ratio_pmx, 1 + ratio_pmx, (ratio_pmx * 2) / rows_needed)]
        df_pmx = pd.concat(pmx_frames, axis=0, ignore_index=True)
        df_pmx = df_pmx.groupby(df_pmx.index // 60).sum()
        # self.save_df_to_plot(df_pmx, 'pmx_simulated_oversampled')

        df = pd.concat([df_nox, df_pmx], axis=1)

        df_combined = pd.concat([df_air.reset_index(), df], axis=1).set_index('time')
        df_combined = df_combined[df_combined.index.notnull()]
        
        # .fillna(method='ffill')
        return df_combined.interpolate(method='time')



    #####################################################################################################
    ################################## DATA COLLECTION FUNCTIONS ########################################
    #####################################################################################################
    ## NOTE: This in crud probably
    async def fetch_simulated_emissions(self, box_id, entries=['CO', 'NOx', 'PMx', 'fuel']):
        # NOTE: First we fetch the simulated emissions
        parser = Parser(self.db, self.sim_id)
        lat, lng = [round(bremicker_boxes[box_id]['lat'], 3), round(bremicker_boxes[box_id]['lng'], 3)]
        raw_emissions = await get_raw_emissions_from_sim(self.db, self.sim_id)
        if raw_emissions == None:
            raw_emissions = parser.parse_emissions()
        # raw_emissions = await parser.get_caqi_data()
        
        # raw_emissions = raw_emissions.reset_index().to_json(orient='index')
        # raw_emissions = raw_emissions.apply(lambda x: x.to_json(orient='records'))
        # raw_emissions = json.loads(raw_emissions)
        # raw_emissions = pd.DataFrame.from_dict(raw_emissions).T
        
        df = pd.DataFrame(pd.read_json(raw_emissions["emissions"], orient='index'))
        # print(df)
        # return raw_emissions
        
        # print(df)
        print(lat, lng)
        # latlng = ['lat', 'lng']
        # df = df[latlng + self.raw_emission_columns]
        # print(df)
        # df = df.groupby([df.index // 60] + latlng)[self.raw_emission_columns].mean().reset_index(latlng)
        # test = df[['time', 'PMx']].set_index('time')
        # print(df.groupby(['time', 'lat', 'lng'])[self.raw_emission_columns].sum())
        # self.save_df_to_plot(test, 'pmx_by_time_simulation')
        mask = (round(df['lat'], 3) == lat) & (round(df['lng'], 3) == lng)
        df = df.loc[mask]
        # print(df)
        
        # df.plot.hist(x='time', y='PMx', grid=True, bins=150, rwidth=0.9, figsize=(12, 8), zorder=2, color='#86bf91')
        # plt.savefig(PLOT_BASEDIR + '/' + 'hist_pmx')
        return df.fillna(method='ffill')
        # return df
        # print(df)
        # plt.hist(df.time, df.NOx)
        
        # df.plot.hist(x='time', y='NOx', grid=True)
        
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


        

    async def fetch_air_and_traffic(self, boxID, start_date='2019-08-01', end_date='2019-11-01', start_hour='7:00', end_hour='10:00'):
        # NOTE: Get real weather data and format it accordingly. Here we'll look at 2019 from 7:00 to 10:00
        df_traffic = await get_bremicker_by_time(
            self.db,
            boxID,
            start_date,
            end_date, 
            start_hour,
            end_hour
        )
        # traffic_mean = df_traffic[boxID].mean()
        # df_traffic = df_traffic.fillna(round(traffic_mean, 2))
        # print(df_traffic)
        df_hawa = await get_hawa_dawa_by_time(
            self.db, 
            start_date,
            end_date, 
            start_hour, 
            end_hour
        )  
        # print(df_hawa)

        df_wind = await fetch_weather_data(start_date, end_date, start_hour, end_hour)
        # print(df_wind)
        df = pd.concat([df_traffic, df_hawa, df_wind], axis=1)
        # print(df)
        # df = df.fillna(method='ffill')
        # return df.fillna(method='bfill')
        return df


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
        # print(df_sum)
        # print(df_sum.mean())
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




    # async def get_air_sensor_data(self):
        # result = await get_all_hawa_dawa(self.db)
        # result = await get_hawa_dawa_by_time(self.db)
        # print(result)
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