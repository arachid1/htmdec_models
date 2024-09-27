from openmsimodel.interactive.gemd_modeller import GEMDModeller
from gemd import MaterialTemplate, ProcessTemplate, MeasurementTemplate, ParameterTemplate, ConditionTemplate, MaterialRun, MaterialSpec, RealBounds, CategoricalBounds, Parameter, NominalReal
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
    print(ni_data['Analysis Parameters'])

    ni_folder_id = find_folder_by_path('66310ca21b56c53abc7e0f9f', ni_data['targetPath'])

    dms_ni_folder_items = client.get(f"/item", parameters={
        'folderId': ni_folder_id,
    })

    area_cag_files = [item['name'] for item in dms_ni_folder_items if item['name'].endswith('_area.cag')]

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
    def make_constant_sr_nanoindentation_experiment(name, strain_rate, ingredient_name):
        nanoindentation_sample_ingredient = Ingredient(ingredient_name)

        tip_area_coefficient_template = ParameterTemplate(
            name="Tip Area Function Coefficient", 
            bounds=RealBounds(0, 100, "")  
        )
        strain_rate_template = ParameterTemplate(
            name="Strain Rate", 
            bounds=RealBounds(0, 100, "mm/s")
        )

        actuator_calibration_coefficient_template = ParameterTemplate(
            name="Actuator Calibration Coefficient", 
            bounds=RealBounds(0, 100, "")
        )

        actuator_column_mass_template = ParameterTemplate(
            name="Actuator Column Mass", 
            bounds=RealBounds(0, 5, "g")
        )

        actuator_spring_stiffness_template = ParameterTemplate(
            name="Actuator Spring Stiffness", 
            bounds=RealBounds(0, 1000, "N/m")
        )

        frame_stiffness_template = ParameterTemplate(
            name="Frame Stiffness", 
            bounds=RealBounds(0, 5000000, "N/m")
        )

        piezo_load_cell_stiffness_template = ParameterTemplate(
            name="Piezo Load Cell Stiffness", 
            bounds=RealBounds(0, 500000, "N/m")
        )

        piezo_calibration_coefficient_template = ParameterTemplate(
            name="Piezo Calibration Coefficient", 
            bounds=RealBounds(0, 100, "")
        )

        piezo_load_cell_damping_coefficient_template = ParameterTemplate(
            name="Piezo Load Cell Damping Coefficient", 
            bounds=RealBounds(0, 100, "")
        )

        piezo_load_cell_effective_mass_template = ParameterTemplate(
            name="Piezo Load Cell Effective Mass", 
            bounds=RealBounds(0, 5, "g")
        )


        nanoindentation_process = Process(
            f"{name} 5-Indent Experiment", 
            template=ProcessTemplate(
                "Nanoindentation Process", 
                parameters=[
                    strain_rate_template,
                    actuator_calibration_coefficient_template,
                    actuator_column_mass_template,
                    actuator_spring_stiffness_template,
                    frame_stiffness_template,
                    piezo_load_cell_stiffness_template,
                    piezo_calibration_coefficient_template,
                    piezo_load_cell_damping_coefficient_template,
                    piezo_load_cell_effective_mass_template,
                    tip_area_coefficient_template
                ]
            )
        )

        nanoindentation_process.update_parameters((Parameter('Strain Rate', value=NominalReal(strain_rate, 'mm/s'))) , which='both')
        for parameter_name in ni_data['Analysis Parameters'].keys():
            unit = 'dimensionless'

            if 'Coefficient' in parameter_name:
                if parameter_name == 'Tip Area Function Coefficients':
                    print(ni_data['Analysis Parameters'][parameter_name])
                    for coefficient_name in ni_data['Analysis Parameters'][parameter_name],keys():
                        percentage = ni_data['Analysis Parameters'][parameter_name][coefficient_name]
                        nanoindentation_process.update_parameters((Parameter(parameter_name, NominalReal(percentage, ''), template=tip_area_coefficient_template)) , which='both')
                elif parameter_name == 'Actuator Calibration Coefficient':
                    template = actuator_calibration_coefficient_template
                elif parameter_name == 'Piezo Calibration Coefficient':
                    template = piezo_calibration_coefficient_template
                elif parameter_name == 'Piezo Load Cell Damping Coefficient':
                    template = piezo_load_cell_damping_coefficient_template

            elif 'Stiffness' in parameter_name:
                unit = 'N/m'
                if parameter_name == 'Actuator Spring Stiffness':
                    template = actuator_spring_stiffness_template
                elif parameter_name == 'Frame Stiffness':
                    template = frame_stiffness_template
                elif parameter_name == 'Pieze Load Cell Stiffness':
                    template = piezo_load_cell_stiffness_template
            elif 'Mass' in parameter_name:
                unit = 'g'
                if parameter_name == 'Actuator Column Mass':
                    template = actuator_column_mass_template
                elif parameter_name == 'Piezo Load Cell Effective Mass':
                    template = piezo_load_cell_effective_mass_template
            param = Parameter(parameter_name, value=NominalReal(ni_data['Analysis Parameters'][parameter_name], unit), template=template)
            # value = NominalReal(ni_data['Analysis Parameters'][parameter_name], unit)
            print("here")
            print(parameter_name)
            print(ni_data['Analysis Parameters'][parameter_name])
            print(type(ni_data['Analysis Parameters'][parameter_name]))
            print(unit)
            nanoindentation_process.update_parameters((param) , which='both')

        indented_sample_material = Material(f"{name} Indented Sample", template=sample_material_template)

        depth_measurement = Measurement(
            name=f"{name} Displacement Measurement", 
            template=MeasurementTemplate(
                name="Depth", 
                parameters=ParameterTemplate(
                    name="Depth", 
                    bounds=RealBounds(0, 10000, "nm")
                )
            )
        )

        force_measurement = Measurement(
            name=f"{name}  Load Measurement", 
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
            name=f"{name} Nanoindentation Experiment on Bonded Sample",
            science_kit=science_kit,  # Reference to the overall experimental kit you're using
            ingredients=[nanoindentation_sample_ingredient],
            material=indented_sample_material,  # The material being tested (attached sample)
            process=nanoindentation_process,  # The nanoindentation process
            measurements=[depth_measurement, force_measurement]  # Measure depth and force during the process
        )

        nanoindentation_sequence.link_within()

        return nanoindentation_sequence

    name = '100'
    ingredient_name = f"{name} SR Attached Sample"
    nanoindentation_calibration_sequence = make_constant_sr_nanoindentation_experiment(f'{name} SR Calibration NI', strain_rate=100, ingredient_name=ingredient_name)
    nanoindentation_calibration_sequence.link_prior(melting_sequence, ingredient_name_to_link=ingredient_name)
    name = '2'
    ingredient_name = f"{name} SR Attached Sample"
    nanoindentation_1_sequence = make_constant_sr_nanoindentation_experiment(f'{name} SR NI', strain_rate=2, ingredient_name=ingredient_name)
    nanoindentation_1_sequence.link_prior(nanoindentation_calibration_sequence, ingredient_name_to_link=ingredient_name)
    name = '1'
    ingredient_name = f"{name} SR Attached Sample"
    nanoindentation_2_sequence = make_constant_sr_nanoindentation_experiment(f'{name} SR NI', strain_rate=1, ingredient_name=ingredient_name)
    nanoindentation_2_sequence.link_prior(nanoindentation_1_sequence, ingredient_name_to_link=ingredient_name)



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
