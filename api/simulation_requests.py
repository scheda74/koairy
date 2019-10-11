import asyncio
from quart import Quart, request, url_for, jsonify
# from app.simulation import parse_emission as parser
from simulation import parse_emission as parser

app = Quart(__name__)

@app.route('/get/emissions')
async def get():
    # TODO: implement method; query database for historic data
    """
    parses files and returns simulated emissions
    """
    emission_cycle = await parser.parse_simulated_emissions()
    return emission_cycle

if __name__ == "__main__":
    app.run('localhost', port=5000, debug=True)