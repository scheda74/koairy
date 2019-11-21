
from app.tools.predictor.neural_nets.py import NeuralNet
from app.tools.predictor.lin_reg import LinReg

from app.db.mongodb import AsyncIOMotorClient, get_database
from fastapi import Depends

import abc


class Predictor(object):
    def __init__(self, context, sim_id=None):
        # self.db = db
        self.sim_id = sim_id
        self.context = context
    
    def predict_emissions(self):
        if self.context == 'linear-regression':
            return LinearRegressionStrategy(self.sim_id).predict_emissions()
        elif self.context == 'lstm':
            return LongShortTermMemoryRecurrentNeuralNetworkStrategy(self.sim_id).predict_emissions()
        elif self.context == 'mlp':
            return MLPRegressorStrategy(self.sim_id).predict_emissions()
        elif self.context == 'cnn':
            print('cnn not yet specified')
            return LinearRegressionStrategy(self.sim_id).predict_emissions()
        else:
            print('Specified strategy not found!')

class PredictorStrategyAbstract(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self, db: AsyncIOMotorClient=Depends(get_database), sim_id=None):
        self.db = db
        self.sim_id = sim_id

    @abc.abstractmethod
    def predict_emissions(self):
        """required method"""

class LinearRegressionStrategy(PredictorStrategyAbstract):
    def predict_emissions(self):
        """Start Linear Regression Model Training and Prediction"""
        lr = LinReg(self.db, self.sim_id)
        return await lr.start_lin_reg()


class MLPRegressorStrategy(PredictorStrategyAbstract):
    def predict_emissions(self):
        """Start MLP Model Training and Prediction"""
        lr = LinReg(self.db, self.sim_id)
        return await lr.start_mlp()

class LongShortTermMemoryRecurrentNeuralNetworkStrategy(PredictorStrategyAbstract):
    def predict_emissions(self):
        """Start LSTM Model Training and Prediction"""
        nn = NeuralNet(self.db, self.sim_id)
        return await nn.start_single_lstm()

class ConvolutionalNeuralNetworkStrategy(PredictorStrategyAbstract):
    def predict_emissions(self):
        """check road and do sth"""
        return {}