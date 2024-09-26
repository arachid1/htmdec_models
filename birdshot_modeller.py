from openmsimodel.interactive.gemd_modeller import GEMDModeller
from gemd import MaterialTemplate, ProcessTemplate, MeasurementTemplate, ParameterTemplate, ConditionTemplate, MaterialRun, MaterialSpec, RealBounds, CategoricalBounds
from openmsimodel.structures.materials_sequence import MaterialsSequence
from openmsimodel.science_kit.science_kit import ScienceKit
from openmsimodel.science_kit.science_kit import ScienceKit
from openmsimodel.entity.gemd.material import Material
from openmsimodel.entity.gemd.process import Process
from openmsimodel.entity.gemd.measurement import Measurement
from openmsimodel.entity.gemd.ingredient import Ingredient
from openmsimodel.structures.materials_sequence import MaterialsSequence
from openmsimodel.db.open_db import OpenDB
from openmsimodel.graph.open_graph import OpenGraph
from openmsimodel.graph.helpers import launch_graph_widget
import re
from girder_client import GirderClient

client = GirderClient(apiUrl='https://data.htmdec.org/api/v1')
client.authenticate(apiKey='MFfpVN81hmOaUV7cTGsovnzdr0iB87ygR0RxkDYA')

def ni_model(file_name, file_path, component):

    print(file_name)
    print(file_path)
    file_id = component['file_id_regex_pattern']
    # print(file_id)
    match = re.search(file_id[0], file_path)

    if match:
        print(f"Matched component: {match.group()}")
    else:
        print("No matched component found")
        return

    raw_data = client.get(
        'entry/search', parameters={'query': f'^{match.group()[:3]}.._VAM-.', 'limit': 1000}
    )
    for ele in raw_data:
        if 'NI-HSR' in ele['data']['targetPath'] and match.group() in ele['data']['targetPath']:
            source = ele['data']
    print(source)

# def match_ni_hsr_folder():
#     match = re.search(file_id[0], file_path)

class BIRDSHOTModeller(GEMDModeller):

    def __init__(self, files_folder, gemd_folder):
        """
        Initialize the GEMDModeller with stores_config, files_folder, and gemd_folder.
        """
        super().__init__(files_folder, gemd_folder)
        self.add_automatable_component(
            lambda s: 'NI-HSR' in s,  
            (r'\b[A-Z]{3}[0-9]{2}\b', True),
            [], 
            lambda file_name, file_path, component: ni_model(file_name, file_path, component)
        )

def main(args=None):
    """
    Main method to run from command line
    """
    BIRDSHOTModeller.run_from_command_line(args)

if __name__ == "__main__":
    main()
