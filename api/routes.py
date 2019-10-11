import asyncio
from quart import Quart, request, url_for, jsonify
# from app.simulation import parse_emission as parser
from simulation import parse_emission as parser
from simulation import simulate_emission as simulator

app = Quart(__name__)

@app.route('/get/emissions')
async def get_emissions():
    # TODO: implement method; query database for historic data
    """
    parses files and returns simulated emissions
    """
    # emission_cycle = await parser.parse_simulated_emissions()
    emission_cycle = await parser.parse_emissions()
    return emission_cycle

@app.route('/get/aqi')
async def get_aqi():
    # TODO: implement method; query database for historic data
    """
    parses files and returns simulated emissions
    """
    # emission_cycle = await parser.parse_simulated_emissions()
    return await parser.get_caqi_data()

@app.route('/start/simulation')
async def start_simulation():
    # TODO: implement method; query database for historic data
    """
    parses files and returns simulated emissions
    """
    # simulator.start()
    print("now I'm parsing results")
    # return await parser.parse_simulated_emissions()
    return await parser.parse_emissions()
    # return "OK"

if __name__ == "__main__":
    app.run('localhost', port=5000, debug=True)