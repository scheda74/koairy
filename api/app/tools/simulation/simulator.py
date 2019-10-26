
import os
import sys
import argparse
import uuid
import datetime
import requests
import xmltodict
import subprocess
import json
import asyncio
from lxml import etree
from app.core.config import EMISSION_OUTPUT_BASE, SUMO_COMMANDLINE

# export SUMO_HOME="/usr/local/opt/sumo/share/sumo"

if 'SUMO_HOME' in os.environ:
    tools = os.path.join(os.environ['SUMO_HOME'], 'tools')
    sys.path.append(tools)
else:
    # os.system("export SUMO_HOME='/usr/local/opt/sumo/share/sumo'")
    # subprocess.call("export SUMO_HOME='/usr/local/opt/sumo/share/sumo'".split())
    sys.exit("please declare environment variable 'SUMO_HOME'")

import traci

class Simulator:
    def __init__(self, cfg_filepath, sim_id):
        self.sim_id = sim_id
        self.cfg_filepath = cfg_filepath
        self.tripinfo_filepath = EMISSION_OUTPUT_BASE + 'tripinfo_%s.xml' % self.sim_id
        self.fcdoutput_filepath = EMISSION_OUTPUT_BASE + 'fcdoutput_%s.xml' % self.sim_id
        self.emission_output_filepath = EMISSION_OUTPUT_BASE + "emission_output_%s.xml" % self.sim_id

    def run(self):
        # step = 0
        while traci.simulation.getMinExpectedNumber() > 0:
            traci.simulationStep()
            # print(step)
            # step += 1
        traci.close()
        sys.stdout.flush()

    async def start(self):
        sumoBinary = SUMO_COMMANDLINE
        sumoCMD = [
            sumoBinary, 
            "-c", self.cfg_filepath,
            # "-c", self.cfg_filepath,
            "--tripinfo-output", self.tripinfo_filepath,
            '--fcd-output', self.fcdoutput_filepath, 
            "--emission-output", self.emission_output_filepath
        ]
        #  "--amitran-output", args.trippath + "trajoutput.xml", 
        # "--phemlight-path", sumo_root + "/data/emissions/PHEMlight/",
        # "--additional-files", create_xml_file(args.lanepath, args.freq, simulation_id), 
        if not os.path.exists(self.emission_output_filepath):
            print(sumoCMD)
            traci.start(sumoCMD, 4041)
            self.run()
        else:
            print("[SIMULATOR] Same simulation already exists. Parsing old file...")
        # doc = {}
        # with open(args.output + "emission_output.xml") as fd:
        #     doc = xmltodict.parse(fd.read())
        # print(doc)
        # return json.dumps(doc)
        
        return