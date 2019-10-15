import asyncio
from quart import Quart, request, url_for, jsonify
# from app.simulation import parse_emission as parser
from app.simulation import parse_emission as parser
from app.simulation.simulator import Simulator
from app.simulation import preprocessor as ip

app = Quart(__name__)

@app.route('/get/emissions')
async def get_emissions():
    # TODO: implement method; query database for historic data
    """
    parses files and returns simulated emissions
    """
    # emission_cycle = await parser.parse_simulated_emissions()
    # emission_cycle = await parser.parse_emissions()
    # return emission_cycle
    return 'Emission'

@app.route('/get/caqi')
async def get_aqi():
    """
    parses files and returns simulated emissions
    """
    # emission_cycle = await parser.parse_simulated_emissions()
    return parser.get_caqi_data()

@app.route('/generate/weights')
async def generate_weights():
    """
    parses files and returns simulated emissions
    """
    processor = ip.PreProcessor()
    await processor.write_weight_file()
    return "File written"

@app.route('/start/simulation')
async def start_simulation():
    """
    parses files and returns simulated emissions
    """
    print("Starting PreProcessor...")
    processor = ip.PreProcessor()
    cfg_filepath = await processor.preprocess_simulation_input()

    print("Starting SUMO...")
    simulator = Simulator(cfg_filepath)
    await simulator.start()
    
    print("Parsing results...")
    return parser.get_caqi_data()
    # return await parser.parse_emissions()
    # return "OK"

if __name__ == "__main__":
    app.run('localhost', port=5000, debug=True)