from typing import List, Optional
from datetime import datetime
import pandas as pd
import requests
import datetime
import calendar
import json
from json import JSONDecodeError
from bson.json_util import dumps
from dateutil.relativedelta import relativedelta

from app.db.mongodb import AsyncIOMotorClient

from app.core.config import (
    database_name,
    BREMICKER_API_KEY,
    BREMICKER_URL,
    bremicker_collection_name
)

async def get_bremicker(conn: AsyncIOMotorClient, start_date='2019-09-01', end_date='2019-09-30'):
    start_date = datetime.datetime.strptime(start_date, '%Y-%m-%d').date()
    end_date = datetime.datetime.strptime(end_date, '%Y-%m-%d').date()
    data = await conn[database_name][bremicker_collection_name].find_one({}, projection={"_id": False})
    if not data:
        await fetch_bremicker(conn)
        # await get_bremicker(conn, start_date.strftime('%Y-%m-%d'), end_date.strftime('%Y-%m-%d'))
    print(data)

    # return await format_to_df(result)



### NOTE: Fetch data from bremicker
async def fetch_bremicker(conn: AsyncIOMotorClient, start_date='2019-10-29', end_date=None):
    print('[MongoDB] No bremicker data found. Fetching from server...')
    start_date = datetime.datetime.strptime(start_date, '%Y-%m-%d').date()
    if end_date:
        end_date = datetime.datetime.strptime(end_date, '%Y-%m-%d').date()
    else:
        end_date = datetime.date.today()
    params = {
        'key': BREMICKER_API_KEY,
        'from': start_date, 
        'to': end_date, 
    }
    try:
        response = requests.get(BREMICKER_URL, params=params)
        response = response.json()
        print('response')
        # print(response)
        await format_bremicker(response)
        # await insert_bremicker(conn, start_date, response)
        
    except JSONDecodeError as error:
        print('error in json decoding: %s' % error)
    else:
        print('all good')
        return response
    
async def format_bremicker(data):
    df = pd.DataFrame(data)
    df['time'] = df[['date', 'time']].apply(lambda x: pd.to_datetime(' '.join(x)), axis=1)
    df = df[['time', 'boxID', 'averageVelocity', 'entryVelocity']]
    df = df.groupby(['boxID', pd.Grouper(key='time', freq='H')]).size()
    # times = df['time']
    # df = df.groupby(['boxID', times.dt.date, times.dt.hour])[['averageVelocity', 'entryVelocity']].mean()
    # [['averageVelocity', 'entryVelocity']].mean()
    # df['veh'] = df.resample('H').transform('count')
    # df = df.groupby(by=df['time'].dt.hour)[['averageVelocity', 'entryVelocity']].mean()
    # df.groupby(['boxID'])[['averageVelocity', 'entryVelocity']]
    # df = df[['date', 'time']]
    print(df)
    # for measure in data:
    #     time = pd.to_datetime(measure['date'] + measure['time'])

async def insert_bremicker(conn: AsyncIOMotorClient, date, data: dict):
    print("[MONGODB] Saving Bremicker data")
    raw_doc = {}
    raw_doc["measure_date"] = date.strftime('%Y-%m-%d')
    raw_doc["data"] = data
    await conn[database_name][bremicker_collection_name].insert_one(raw_doc)
    print("[MONGODB] Bremicker data successfully saved!")