# import pymongo

# class DB(object):
#     URI = "mongodb://localhost:27017"

#     @staticmethod
#     def init():
#         client = pymongo.MongoClient(DB.URI)
#         DB.DATABASE = client['mongo-db']

#     @staticmethod
#     def insert(collection, data):
#         DB.DATABASE[collection].insert(data)

#     @staticmethod
#     def find_one(collection, query):
#         return DB.DATABASE[collection].find_one(query)


from motor.motor_asyncio import AsyncIOMotorClient


class DataBase:
    client: AsyncIOMotorClient = None


db = DataBase()


async def get_database() -> AsyncIOMotorClient:
    return db.client