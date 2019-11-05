


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
        return None
        # raise RuntimeError(f" Couldn't find caqi emissions for specified simulation,"
        #                    f" sim_id={sim_id} emission_id={emission_doc}")

async def get_raw_emissions_from_sim(conn: AsyncIOMotorClient, sim_id: str):
    emission_doc = await conn[database_name][raw_emission_collection_name].find_one({"sim_id": sim_id}, projection={"id": False})
    if emission_doc:
        return emission_doc
    else:
        return None
        # raise RuntimeError(f" Couldn't find raw emissions data for specified simulation,"
        #                    f" sim_id={sim_id} emission_id={emission_doc}")


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