import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import datetime
import json
import math

from fastapi import Depends
from app.db.mongodb import AsyncIOMotorClient, get_database

from keras.models import Sequential
from keras.layers import Dense
from keras.layers import LSTM
from sklearn.linear_model import LinearRegression
from sklearn.metrics import (mean_absolute_error, mean_squared_error)
from sklearn.neural_network import MLPRegressor
from sklearn import preprocessing as pre

from app.tools.predictor.utils.model_preprocessor import ModelPreProcessor

class NeuralNet():
    def __init__(self, db: AsyncIOMotorClient, sim_id=None, existing_regr=None):
        self.db = db
        self.sim_id = sim_id
        self.existing_regr = existing_regr
        self.raw_emission_columns = ['CO', 'NOx', 'PMx']
        self.real_emission_columns = ['no2', 'pm2.5', 'pm10', 'o3', 'WIND_SPEED', 'WIND_DIR']
        self.mp = ModelPreProcessor(self.db, self.sim_id)

    async def start_single_lstm(
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
        df = await self.mp.aggregate_data(boxID, start_date, end_date, start_hour, end_hour)
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
        trainX, trainY = self.mp.create_dataset(df_train[[output_key]].values, look_back)
        testX, testY = self.mp.create_dataset(df_test[[output_key]].values, look_back)

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
        print("")
