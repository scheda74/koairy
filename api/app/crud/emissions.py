


from typing import List, Optional
from datetime import datetime
import requests
import datetime
import calendar
import json
from bson.json_util import dumps
from dateutil.relativedelta import relativedelta

from app.db.mongodb import AsyncIOMotorClient
from app.core.config import (
    database_name, 
    caqi_emission_collection_name,
    raw_emission_collection_name,
    bremicker_collection_name,
    air_hawa_collection_name,
    HAWA_DAWA_URL,
    HAWA_DAWA_API_KEY
)


async def get_caqi_emissions_for_sim(conn: AsyncIOMotorClient, sim_id: str):
    emission_doc = await conn[database_name][caqi_emission_collection_name].find_one({"sim_id": sim_id}, projection={"id": False})
    if emission_doc:
        return emission_doc
    else:
        raise RuntimeError(f" Couldn't find caqi emissions for specified simulation,"
                           f" sim_id={sim_id} emission_id={emission_doc}")

async def get_raw_emissions_from_sim(conn: AsyncIOMotorClient, sim_id: str):
    emission_doc = await conn[database_name][raw_emission_collection_name].find_one({"sim_id": sim_id}, projection={"id": False})
    if emission_doc:
        return emission_doc
    else:
        raise RuntimeError(f" Couldn't find raw emissions data for specified simulation,"
                           f" sim_id={sim_id} emission_id={emission_doc}")


async def insert_caqi_emissions(conn: AsyncIOMotorClient, sim_id: str, emissions: dict):
    print("[MONGODB] Saving calculated CAQI emissions")
    caqi_doc = {}
    caqi_doc["emissions"] = emissions
    caqi_doc["created_at"] = datetime.datetime.utcnow()
    caqi_doc["sim_id"] = sim_id
    await conn[database_name][caqi_emission_collection_name].insert_one(caqi_doc)

async def insert_raw_emissions(conn: AsyncIOMotorClient, sim_id: str, emissions: dict):
    print("[MONGODB] Saving simulated emissions")
    raw_doc = {}
    raw_doc["emissions"] = emissions
    raw_doc["created_at"] = datetime.datetime.utcnow()
    raw_doc["sim_id"] = sim_id
    await conn[database_name][raw_emission_collection_name].insert_one(raw_doc)

async def insert_air_traffic(conn: AsyncIOMotorClient, date, data: dict):
    print("[MONGODB] Saving HawaDawas data of {date}")
    raw_doc = {}
    raw_doc["measure_date"] = date.strftime('%Y-%m-%d')
    raw_doc["data"] = json.dumps(data)
    await conn[database_name][air_hawa_collection_name].insert_one(raw_doc)

async def get_all_air_traffic(conn: AsyncIOMotorClient):
    result = []
    async for data in conn[database_name][air_hawa_collection_name].find({}, projection={"id": False}):
        if data:
            result.append(data)
            break
        else:
            break    
    
    if len(result) == 0: 
        await fetch_air_traffic(conn, '2019-01-01')
    return result
    
        

async def fetch_air_traffic(conn: AsyncIOMotorClient, date):
    date = datetime.datetime.strptime(date, '%Y-%m-%d').date()
    current_date = datetime.date.today()
    result = []
    while date < current_date:
        first_day, last_day = get_month_day_range(date)
        data = await fetch_air_traffic_from_hawa_dawa(first_day, last_day)
        await insert_air_traffic(conn, first_day, data)
        result.append(data)
        date = date + relativedelta(months=1)
    # print(result)
    return result

async def fetch_air_traffic_from_hawa_dawa(start_date, end_date):
    params = {
        'apikey': HAWA_DAWA_API_KEY,
        'timeseries_startdate': start_date, 
        'timeseries_enddate': end_date, 
        'show_plausible': False,
        'format': 'geojson',
        'missing_value': -999,
        'crs': 'global'
    }
    response = requests.get(HAWA_DAWA_URL, params=params)
    return response.json()    

def get_month_day_range(date):
    first_day = date.replace(day = 1)
    last_day = date.replace(day = calendar.monthrange(date.year, date.month)[1])
    return first_day, last_day