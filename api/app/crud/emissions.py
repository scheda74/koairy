


from typing import List, Optional
from datetime import datetime

from ..db.mongodb import AsyncIOMotorClient
from ..core.config import (
    database_name, 
    caqi_emission_collection_name,
    raw_emission_collection_name
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

# DB.insert(collection='caqi_emissions', data={ "created_at:": datetime.datetime.utcnow(), "sim_id": self.sim_id, "emissions": result })
async def insert_caqi_emissions(conn: AsyncIOMotorClient, sim_id: str, emissions: dict):
    caqi_doc = {}
    caqi_doc["emissions"] = emissions
    caqi_doc["created_at"] = datetime.utcnow()
    caqi_doc["sim_id"] = sim_id
    await conn[database_name][caqi_emission_collection_name].insert_one(caqi_doc)

async def insert_raw_emissions(conn: AsyncIOMotorClient, sim_id: str, emissions: dict):
    raw_doc = {}
    raw_doc["emissions"] = emissions
    raw_doc["created_at"] = datetime.utcnow()
    raw_doc["sim_id"] = sim_id
    await conn[database_name][raw_emission_collection_name].insert_one(raw_doc)