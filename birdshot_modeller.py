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

# Function to get subfolders of a folder using the Girder API
def get_subfolders(parent_folder_id):
    
    # API request to get subfolders of a folder
    subfolders = client.get(f"/folder", parameters={
        'parentType': 'folder',
        'parentId': parent_folder_id
    })
    
    return subfolders

# Function to find a subfolder by its name within a parent folder
def find_subfolder_by_name(parent_folder_id, folder_name):
    subfolders = get_subfolders(parent_folder_id)
    
    # Search for the folder by name
    for folder in subfolders:
        if folder['name'] == folder_name:
            return folder  # Return the folder metadata
    
    return None  # Return None if not found

# Function to find a folder by traversing a given path
def find_folder_by_path(root_folder_id, folder_path):
    folder_names = folder_path.strip("/").split("/")  # Split path into folder names
    
    current_folder_id = root_folder_id
    
    # Traverse through the folder hierarchy based on the folder names
    for folder_name in folder_names:
        folder = find_subfolder_by_name(current_folder_id, folder_name)
        if folder:
            current_folder_id = folder['_id']  # Update current folder ID to the found folder
        else:
            print(f"Folder '{folder_name}' not found.")
            return None
    
    return current_folder_id  # Return the ID of the final folder if found

def ni_model(file_name, file_path, component):

    print(file_name)
    print(file_path)
    file_id = component['file_id_regex_pattern']
    match = re.search(file_id[0], file_path)

    if not match:
        print("No pattern found.")
        return

    raw_data_forms = client.get(
        'entry/search', parameters={'query': f'^{match.group()[:3]}.._VAM-.', 'limit': 1000}
    )
    for form in raw_data_forms:
        if 'NI-HSR' in form['data']['targetPath'] and match.group() in form['data']['targetPath']:
            ni_data = form['data']

    print(ni_data)

    ni_folder_id = find_folder_by_path('66310ca21b56c53abc7e0f9f', ni_data['targetPath'])

    print('folder_id', ni_folder_id)
    dms_ni_folder_items = client.get(f"/item", parameters={
        'folderId': ni_folder_id,
    })
    print(dms_ni_folder_items)

    area_cag_files = [item['name'] for item in dms_ni_folder_items if item['name'].endswith('_area.cag')]
    print(area_cag_files)

    science_kit = ScienceKit()

    ######
    sample_material_template = MaterialTemplate("Sample")
    temperature_measurement_template = MeasurementTemplate("Temperature")
    crystal_bond_ingredient = Ingredient("Crystal Bond")
    sample_ingredient = Ingredient("Sample")

    melting_process = Process(
        "Melting Crystal Bond", 
        template=ProcessTemplate(
            "Melting", 
            parameters=ParameterTemplate(
                name="Temperature",
                bounds=RealBounds(50, 200, "Celsius"),  # Assuming melting happens between 50°C to 200°C
            )
        )
    )

    # Define the resulting material (sample attached with crystal bond)
    attached_sample_material = Material("Sample Attached with Crystal Bond", template=sample_material_template)

    # Define a measurement, for example, the temperature during melting
    melting_temperature_measurement = Measurement(
        'Melting Temperature', 
        template=temperature_measurement_template
    )

    # Define the melting sequence using MaterialsSequence
    melting_sequence = MaterialsSequence(
        name="Melting Crystal Bond to Attach Sample",
        science_kit=science_kit,
        material=attached_sample_material,
        ingredients=[crystal_bond_ingredient, sample_ingredient],
        process=melting_process,
        measurements=[melting_temperature_measurement]
    )

    # Link internal elements within the sequence
    melting_sequence.link_within()

    #####
    nanoindentation_sample_ingredient = Ingredient("Sample Attached with Crystal Bond")

    nanoindentation_process = Process(
        "Nanoindentation", 
        template=ProcessTemplate(
            "Nanoindentation Process", 
            parameters=[
                ParameterTemplate(
                    name="Indentation Depth", 
                    bounds=RealBounds(0, 10000, "nm")  # Depth range from 0 to 10,000 nm
                ),
                ParameterTemplate(
                    name="Applied Force", 
                    bounds=RealBounds(0, 0.5, "micronewton")  # Force applied in micronewtons
                )
            ]
        )
    )

    indented_sample_material = Material("Indented Sample", template=sample_material_template)

    depth_measurement = Measurement(
        name="Indentation Depth Measurement", 
        template=MeasurementTemplate(
            name="Depth", 
            parameters=ParameterTemplate(
                name="Depth", 
                bounds=RealBounds(0, 10000, "nm")
            )
        )
    )

    force_measurement = Measurement(
        name="Applied Force Measurement", 
        template=MeasurementTemplate(
            name="Force", 
            parameters=ParameterTemplate(
                name="Force", 
                bounds=RealBounds(0, 500, "micronewton")
            )
        )
    )

    # Define the nanoindentation sequence, linking it to the sample and adding measurements
    nanoindentation_sequence = MaterialsSequence(
        name="Nanoindentation Experiment on Bonded Sample",
        science_kit=science_kit,  # Reference to the overall experimental kit you're using
        ingredients=[nanoindentation_sample_ingredient],
        material=indented_sample_material,  # The material being tested (attached sample)
        process=nanoindentation_process,  # The nanoindentation process
        measurements=[depth_measurement, force_measurement]  # Measure depth and force during the process
    )

    nanoindentation_sequence.link_within()

    nanoindentation_sequence.link_prior(melting_sequence, ingredient_name_to_link="Sample Attached with Crystal Bond")

    return science_kit.assets()








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
