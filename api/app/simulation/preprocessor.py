import os, sys
import subprocess
import xml.etree.ElementTree as ET
from lxml import etree
from app.simulation.constants import Constants


# path_to_net,
# path_to_trip_output,
# simulation_timesteps,
# path_to_route_output,
# path_to_weights,
# '../data/input-simulation/kirchheim-street-names.net.xml',
# '../data/input-simulation/scenario{0}.trips.xml'.format(scenario_id),
# '../data/input-simulation/scenario{0}.rou.xml'.format(scenario_id),
# '../data/input-simulation/scenario{0}'.format(scenario_id),


class PreProcessor(Constants):
    def __init__(
        self,
        new_net_path=None, 
        timesteps=10800, 
        agents=9500,
        weights_per_area={
            'aschheim_west': 0.1,
            'ebersberg_east': 0.37,
            'feldkirchen_west': 0.1,
            'heimstetten_industrial_1': 0.01,
            'heimstetten_industrial_2': 0.01,
            'heimstetten_residential': 0.18,
            'kirchheim_industrial_east': 0.1,
            'kirchheim_industrial_west': 0.1,
            'kirchheim_residential': 0.18,
            'unassigned_edges': 0.05
        },
        weight_type='src',
        fringe_factor=1
    ):
        self.timesteps = timesteps
        self.agents = agents
        self.weights_per_area = weights_per_area
        self.weight_type = weight_type
        self.fringe_factor = fringe_factor
        self.new_net_path = new_net_path

        weights = "".join([str(v).replace('.', '') for v in self.weights_per_area.values()])
        self.weights_filepath = Constants.WEIGHT_INPUT + "%s%s-%s.%s.xml" % (self.agents, self.timesteps, weights, self.weight_type)
        self.trip_filepath = Constants.TRIP_OUTPUT + weights + ".trip.xml"
        self.route_filepath = Constants.ROUTE_OUTPUT + weights + ".rou.xml"
        self.net_filepath = new_net_path if new_net_path != None else Constants.DEFAULT_NET_INPUT
        self.cfg_filepath = Constants.SUMO_CFG + weights + ".sumocfg"

    def write_sumocfg_file(self):
        configuration_tag = ET.Element('configuration')

        input_tag = ET.SubElement(
            configuration_tag,
            'input'
        )
        net_file_tag = ET.SubElement(
            input_tag,
            'net-file',
            {'value': self.net_filepath}
        )
        route_files_tag = ET.SubElement(
            input_tag,
            'route-files',
            {'value': self.route_filepath}
        )
        gui_settings_file_tag = ET.SubElement(
            input_tag,
            'gui-settings-file',
            {'value': 'gui-settings-origin-dest-vehicles.cfg'}
        )

        time_tag = ET.SubElement(configuration_tag, 'time')
        begin_tag = ET.SubElement(time_tag, 'begin', {'value': '0'})
        end_tag = ET.SubElement(time_tag, 'end', {'value': str(self.timesteps)})

        processing_tag = ET.SubElement(configuration_tag, 'processing')
        junction_blocker_ignore_tag = ET.SubElement(processing_tag, 'ignore-junction-blocker', {'value': '1'})

        configuration_file = ET.ElementTree(configuration_tag)
        configuration_file.write(self.cfg_filepath)

        # pretty formatting
        parser = etree.XMLParser(remove_blank_text=True)
        tree = etree.parse(self.cfg_filepath, parser)
        tree.write(self.cfg_filepath, pretty_print=True)

    def init_weight_file(self, edges_per_area):
        print("Weight file not found - Initializing new weight file...")
        root = ET.Element('edgedata')
        interval = ET.SubElement(root, 'interval', {'begin': '0', 'end': str(self.timesteps)})

        for area_id, edges in edges_per_area.items():
            for edge in edges:
                ET.SubElement(interval, 'edge', {'id': edge, 'value': '0'})

        tree = ET.ElementTree(root)
        tree.write(self.weights_filepath)

        #pretty formatting
        parser = etree.XMLParser(remove_blank_text=True)
        tree = etree.parse(self.weights_filepath, parser)
        tree.write(self.weights_filepath, pretty_print=True)


    def write_weight_file(self):
        print("Writing/Formatting weight file...")
        weights_per_area_keys = set(self.weights_per_area.keys())

        if (weights_per_area_keys.symmetric_difference(Constants.VALID_AREA_IDS) != set()):
            raise ValueError('area_ids must only be exactly %r.' % Constants.VALID_AREA_IDS)
        if (self.weight_type != 'src' and self.weight_type != 'dst'):
            raise ValueError("weight_type must be 'src' or 'dst'")

        edges_per_area = {}
        taz_tree = ET.parse(Constants.AREA_OF_INTEREST)
        taz_root = taz_tree.getroot()

        for taz in taz_root.iter('taz'):
            if (taz.get('id') in Constants.VALID_AREA_IDS):
                edges_per_area['{0}'.format(taz.get('id'))] = taz.get('edges').split()
        
        if not os.path.exists(self.weights_filepath):
            self.init_weight_file(edges_per_area)
        
        scenario_weights_tree = ET.parse(self.weights_filepath)
        scenario_weights_root = scenario_weights_tree.getroot()

        for edge in scenario_weights_root.iter('edge'):
            edge_id = edge.get('id')
            if (edge_id in edges_per_area['aschheim_west']):
                edge.set('value', '{0}'\
                .format(self.weights_per_area['aschheim_west']/len(edges_per_area['aschheim_west'])))
            elif (edge_id in edges_per_area['ebersberg_east']):
                edge.set('value', '{0}'\
                .format(self.weights_per_area['ebersberg_east']/len(edges_per_area['ebersberg_east'])))
            elif (edge_id in edges_per_area['feldkirchen_west']):
                edge.set('value', '{0}'\
                .format(self.weights_per_area['feldkirchen_west']/len(edges_per_area['feldkirchen_west'])))
            elif (edge_id in edges_per_area['heimstetten_industrial_1']):
                edge.set('value', '{0}'\
                .format(self.weights_per_area['heimstetten_industrial_1']/len(edges_per_area['heimstetten_industrial_1'])))
            elif (edge_id in edges_per_area['heimstetten_industrial_2']):
                edge.set('value', '{0}'\
                .format(self.weights_per_area['heimstetten_industrial_2']/len(edges_per_area['heimstetten_industrial_2'])))
            elif (edge_id in edges_per_area['heimstetten_residential']):
                edge.set('value', '{0}'\
                .format(self.weights_per_area['heimstetten_residential']/len(edges_per_area['heimstetten_residential'])))
            elif (edge_id in edges_per_area['kirchheim_industrial_east']):
                edge.set('value', '{0}'\
                .format(self.weights_per_area['kirchheim_industrial_east']/len(edges_per_area['kirchheim_industrial_east'])))
            elif (edge_id in edges_per_area['kirchheim_industrial_west']):
                edge.set('value', '{0}'\
                .format(self.weights_per_area['kirchheim_industrial_west']/len(edges_per_area['kirchheim_industrial_west'])))
            elif (edge_id in edges_per_area['kirchheim_residential']):
                edge.set('value', '{0}'\
                .format(self.weights_per_area['kirchheim_residential']/len(edges_per_area['kirchheim_residential'])))
            elif (edge_id in edges_per_area['unassigned_edges']):
                edge.set('value', '{0}'\
                .format(self.weights_per_area['unassigned_edges']/len(edges_per_area['unassigned_edges'])))

        scenario_weights_tree.write(self.weights_filepath)
        return

    def write_random_trips_and_routes(self):
        # path_to_script = tools + '/randomTrips.py'
        print("Writing random trips and route files...")
        cmd = "python %s -n %s -e %s -o %s --route-file %s --validate --fringe-factor %s -p %s --weights-prefix %s"\
                % ( Constants.RANDOM_TRIP_TOOL,
                    self.net_filepath,
                    self.timesteps, 
                    self.trip_filepath, 
                    self.route_filepath,
                    self.fringe_factor, 
                    ((self.timesteps-0) / (self.agents * 1.0)),
                    self.weights_filepath
                )
        subprocess.call(cmd.split())

    async def preprocess_simulation_input(self):
        # if self.new_net_path is not None:
        #     write_taz_file(
        #         '../data/input-simulation/kirchheim-street-names.net.xml',
        #         '../data/input-simulation/areas-of-interest.boundaries.xml',
        #         '../data/input-simulation/areas-of-interest.taz.xml'
        #     )
        if not os.path.exists(self.weights_filepath):
            self.write_weight_file()
            self.write_random_trips_and_routes()
            self.write_sumocfg_file()
        
        return self.cfg_filepath
