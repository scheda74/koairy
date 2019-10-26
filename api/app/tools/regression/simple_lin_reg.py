import numpy as np 
import pandas as pd
import json

from fastapi import Depends
from app.db.mongodb import AsyncIOMotorClient, get_database
from app.crud.emissions import get_raw_emissions_from_sim


class LinReg():
    def __init__(self, db: AsyncIOMotorClient, sim_id):
        self.db = db
        self.sim_id = sim_id
        self.columns = ['CO', 'NOx', 'PMx']
    
    async def fetch_simulated_emissions(self):
        raw_emissions = await get_raw_emissions_from_sim(self.db, self.sim_id)
        # print(raw_emissions["emissions"])
        pd_emissions = pd.read_json(raw_emissions["emissions"], orient='index')
        # df = pd.DataFrame(pd_emissions, dtype=float)
        df = pd.DataFrame(pd_emissions)
        print(df[self.columns])
        return raw_emissions
