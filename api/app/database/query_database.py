from app.database.database import DB
import json
from bson import ObjectId

def get_latest_emissions():
    collection = DB.DATABASE["caqi_emissions"]
    cursor = list(collection.find({}, {"_id": False}).sort([("created_at", -1)]).limit(1))
    if len(cursor):
        data = cursor[0]
        # data["created_at"] = str(data["created_at"])
        return data
    else:
        return None