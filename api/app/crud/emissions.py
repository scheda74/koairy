


from typing import List, Optional
from bson import ObjectId
from slugify import slugify
from datetime import datetime

from ..db.mongodb import AsyncIOMotorClient
from ..core.config import database_name, emission_collection_name


async def get_caqi_emissions_for_sim(conn: AsyncIOMotorClient, sim_id: str):
    emission_doc = await conn[database_name][emission_collection_name].find_one({"sim_id": sim_id}. projection={"id": False})
    if emission_doc:
        return emission_doc["emissions"]
    else:
        raise RuntimeError(f" Couldn't find caqi emissions for specified simulation,"
                           f" sim_id={sim_id} emission_id={emission_doc}")