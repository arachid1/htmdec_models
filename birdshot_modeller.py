from openmsimodel.interactive.gemd_modeller import GEMDModeller
from gemd import MaterialTemplate, ProcessTemplate, MeasurementTemplate, ParameterTemplate, ConditionTemplate, MaterialRun, MaterialSpec, RealBounds, CategoricalBounds, Parameter, NominalReal, FileLink
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

science_kit = ScienceKit()

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
                    for coefficient_name in ni_data['Analysis Parameters'][parameter_name].keys():
                        percentage = ni_data['Analysis Parameters'][parameter_name][coefficient_name]
                        print("aqui")
                        print(percentage)
                        print(type(percentage))
                        param = Parameter(parameter_name, value=NominalReal(percentage, unit), template=tip_area_coefficient_template)
                        nanoindentation_process.update_parameters((param), which='both')
                    continue
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

            # value = NominalReal(ni_data['Analysis Parameters'][parameter_name], unit)
            print("here")
            print(parameter_name)
            print(ni_data['Analysis Parameters'][parameter_name])
            print(type(ni_data['Analysis Parameters'][parameter_name]))
            print(unit)
            param = Parameter(parameter_name, value=NominalReal(ni_data['Analysis Parameters'][parameter_name], unit), template=template)

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
            name=f"{name} Load Measurement", 
            template=MeasurementTemplate(
                name="Force", 
                parameters=ParameterTemplate(
                    name="Force", 
                    bounds=RealBounds(0, 500, "micronewton")
                )
            )
        )

        hardness_measurement = Measurement(
            name=f"{name} Hardness Measurement", 
            template=MeasurementTemplate(
                name="Hardness", 
                parameters=ParameterTemplate(
                    name="Hardness", 
                    bounds=RealBounds(-500, 500, "GPa")
                )
            )
        )

        area_measurement = Measurement(
            name=f"{name} Area Measurement", 
            template=MeasurementTemplate(
                name="Area", 
                parameters=ParameterTemplate(
                    name="Area", 
                    bounds=RealBounds(-500, 500, "nm^2")
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
            measurements=[depth_measurement, force_measurement, hardness_measurement, area_measurement]  # Measure depth and force during the process
        )

        nanoindentation_sequence.link_within()

        return nanoindentation_sequence

    sr_value = '100'
    ingredient_name = f"Attached Sample ({sr_value} SR)"
    nanoindentation_calibration_sequence = make_constant_sr_nanoindentation_experiment(f'Calibration ({sr_value} SR)', strain_rate=100, ingredient_name=ingredient_name)
    nanoindentation_calibration_sequence.link_prior(melting_sequence, ingredient_name_to_link=ingredient_name)
    csr_files = [(item['name'], item['_id']) for item in dms_ni_folder_items if '_CSR_{}'.format(sr_value) in item['name'] and (not('.vk6' in item['name']))]
    for csr_file in csr_files:
        nanoindentation_calibration_sequence.process.update_filelinks(FileLink(f'{csr_file[0]} (DMS)', url=f'https://data.htmdec.org/api/v1/item/{csr_file[1]}'), which='run')

    sr_value = '2'
    ingredient_name = f"Attached Sample ({sr_value} SR)"
    nanoindentation_1_sequence = make_constant_sr_nanoindentation_experiment(f'({sr_value} SR)', strain_rate=2, ingredient_name=ingredient_name)
    nanoindentation_1_sequence.link_prior(nanoindentation_calibration_sequence, ingredient_name_to_link=ingredient_name)
    csr_files = [(item['name'], item['_id']) for item in dms_ni_folder_items if '_CSR_{}'.format(sr_value) in item['name'] and (not('.vk6' in item['name']))]
    for csr_file in csr_files:
        nanoindentation_1_sequence.process.update_filelinks(FileLink(f'{csr_file[0]} (DMS)', url=f'https://data.htmdec.org/api/v1/item/{csr_file[1]}'), which='run')

    sr_value = '1'
    ingredient_name = f"Attached Sample ({sr_value} SR)"
    nanoindentation_2_sequence = make_constant_sr_nanoindentation_experiment(f'({sr_value} SR)', strain_rate=1, ingredient_name=ingredient_name)
    nanoindentation_2_sequence.link_prior(nanoindentation_1_sequence, ingredient_name_to_link=ingredient_name)
    csr_files = [(item['name'], item['_id']) for item in dms_ni_folder_items if (('_CSR_{}'.format(sr_value) in item['name']) and (not('.vk6' in item['name'])) and (not ('_CSR_100' in item['name']))) ]
    for csr_file in csr_files:
        nanoindentation_2_sequence.process.update_filelinks(FileLink(f'{csr_file[0]} (DMS)', url=f'https://data.htmdec.org/api/v1/item/{csr_file[1]}'), which='run')

    area_cag_files = [item['name'] for item in dms_ni_folder_items if item['name'].endswith('_area.cag')]

    ######
    def make_profilometry_experiment(name, ingredient_name):
        """
        Create a MaterialsSequence for Profilometry using the Keyence Microscope
        to capture images and measurements of the indents.
        """
        profilometry_sample_ingredient = Ingredient(ingredient_name)

        # Define the parameter templates for profilometry
        projected_contact_area_template = ParameterTemplate(
            name="Projected Contact Area",
            bounds=RealBounds(0, 1000, "µm²")
        )
        
        pile_up_ratio_template = ParameterTemplate(
            name="Pile-up Ratio",
            bounds=RealBounds(0, 1, "dimensionless")
        )

        # Define the profilometry process
        profilometry_process = Process(
            f"Keyence Profilometry Step {name}", 
            template=ProcessTemplate(
                "Profilometry Process",
                parameters=[
                    projected_contact_area_template,
                    pile_up_ratio_template
                ]
            )
        )

        # Create a measurement for 3D scans
        scan_measurement = Measurement(
            name=f"3D Indent Scan {name}", 
            template=MeasurementTemplate(
                name="Indent Scan", 
                parameters=ParameterTemplate(
                    name="Indent Depth", 
                    bounds=RealBounds(0, 10000, "nm")  # Measurement of indent depth
                )
            )
        )

        # Define the resulting material (the scanned sample)
        scanned_sample_material = Material(f"Scanned Sample {name}", template=sample_material_template)

        # Define the profilometry sequence using MaterialsSequence
        profilometry_sequence = MaterialsSequence(
            name=f"Profilometry {name} on Indented Sample",
            science_kit=science_kit,
            material=scanned_sample_material,
            ingredients=[profilometry_sample_ingredient],
            process=profilometry_process,
            measurements=[scan_measurement]  # Measure the indent scans
        )

        # Link internal elements within the sequence
        profilometry_sequence.link_within()

        return profilometry_sequence

    # Example usage:
    # Profilometry step after the nanoindentation experiment
    ingredient_name = f"Indented Sample ({sr_value} SR)"
    profilometry_sequence = make_profilometry_experiment(f'1', ingredient_name)
    profilometry_sequence.link_prior(nanoindentation_2_sequence, ingredient_name_to_link=ingredient_name)

    # Now, map the .vk6 file associated with the profilometry scan to the sequence
    vk6_files = [(item['name'], item['_id']) for item in dms_ni_folder_items if item['name'].endswith('.vk6')]

    for vk6_file in vk6_files:
        profilometry_sequence.process.update_filelinks(FileLink(f'{vk6_file[0]} (Profilometry)', url=f'https://data.htmdec.org/api/v1/item/{vk6_file[1]}'), which='run')
    

    return science_kit.assets()


def eds_model(file_name, file_path, component):
    file_id = component['file_id_regex_pattern']
    match = re.search(file_id[0], file_path)

    if not match:
        print("No pattern found.")
        return

    raw_data_forms = client.get(
        'entry/search', parameters={'query': f'^{match.group()[:3]}.._VAM-.', 'limit': 1000}
    )
    for form in raw_data_forms:
        if 'EDS' in form['data']['targetPath'] and match.group() in form['data']['targetPath']:
            ebsd_eds_data = form['data']
    # print(ebsd_eds_data)
    def make_ebsd_eds_mapping_sequence(data):

        # Define material template for the sample
        sample_material_template = MaterialTemplate("Sample")

        # Ingredients (e.g., the material being analyzed)
        sample_ingredient = Ingredient("Sample")

        # Create process template for EBSD and EDS Mapping
        ebsd_eds_mapping_process_template = ProcessTemplate(
            name="EBSD and EDS Mapping",
            parameters=[
                ParameterTemplate(name="Beam Current", bounds=RealBounds(-1, 100, "nanoampere")),
                ParameterTemplate(name="Beam Voltage", bounds=RealBounds(0, 30, "kilovolt")),
                ParameterTemplate(name="Dwell Time", bounds=RealBounds(-1, 10, "seconds")),
                ParameterTemplate(name="Sample Tilt", bounds=RealBounds(-1, 90, "degrees")),
                ParameterTemplate(name="Working Distance", bounds=RealBounds(-1, 20, "millimeter")),
                ParameterTemplate(name="Low Vacuum", bounds=CategoricalBounds(["None", "Low", "Medium", "High"]))
            ]
        )

        # Create EBSD and EDS Mapping process
        ebsd_eds_mapping_process = Process(
            name="EBSD and EDS Mapping",
            template=ebsd_eds_mapping_process_template
        )

        # Fill EBSD and EDS Mapping process with values from the form
        ebsd_eds_mapping_process.update_parameters(
            Parameter('Beam Current', NominalReal(data['EBSD and EDS Mapping']['Beam Current'], 'nanoampere'), template=ebsd_eds_mapping_process_template.parameters[0]),
            Parameter('Beam Voltage', NominalReal(data['EBSD and EDS Mapping']['Beam Voltage'], 'kilovolt'), template=ebsd_eds_mapping_process_template.parameters[1]),
            Parameter('Dwell Time', NominalReal(data['EBSD and EDS Mapping']['Dwell Time'], 'seconds'), template=ebsd_eds_mapping_process_template.parameters[2]),
            Parameter('Sample Tilt', NominalReal(data['EBSD and EDS Mapping']['Sample Tilt'], 'degrees'), template=ebsd_eds_mapping_process_template.parameters[3]),
            Parameter('Working Distance', NominalReal(data['EBSD and EDS Mapping']['Working Distance'], 'millimeter'), template=ebsd_eds_mapping_process_template.parameters[4]),
            Parameter('Low Vacuum', CategoricalValue(data['EBSD and EDS Mapping']['Low Vacuum']), template=ebsd_eds_mapping_process_template.parameters[5])
        )

        # Create Material for the EBSD and EDS Mapping sample
        ebsd_sample_material = Material("Sample from EBSD and EDS Mapping", template=sample_material_template)

        # Create Measurement templates for EDS Measured Composition and StdDev
        eds_composition_measurement_template = MeasurementTemplate("EDS Measured Composition")
        eds_stddev_measurement_template = MeasurementTemplate("EDS Measured Composition StdDev")

        # Create Measurements for EDS composition results (primary phase)
        eds_measured_composition = Measurement(
            name="EDS Measured Composition",
            template=eds_composition_measurement_template,
            parameters=[
                Parameter('Al', NominalReal(data['Results']['Measured Composition (%)']['Al'], '%'), template=eds_composition_measurement_template),
                Parameter('Co', NominalReal(data['Results']['Measured Composition (%)']['Co'], '%'), template=eds_composition_measurement_template),
                Parameter('Cr', NominalReal(data['Results']['Measured Composition (%)']['Cr'], '%'), template=eds_composition_measurement_template),
                Parameter('Cu', NominalReal(data['Results']['Measured Composition (%)']['Cu'], '%'), template=eds_composition_measurement_template),
                Parameter('Fe', NominalReal(data['Results']['Measured Composition (%)']['Fe'], '%'), template=eds_composition_measurement_template),
                Parameter('Mn', NominalReal(data['Results']['Measured Composition (%)']['Mn'], '%'), template=eds_composition_measurement_template),
                Parameter('Ni', NominalReal(data['Results']['Measured Composition (%)']['Ni'], '%'), template=eds_composition_measurement_template),
                Parameter('V', NominalReal(data['Results']['Measured Composition (%)']['V'], '%'), template=eds_composition_measurement_template)
            ]
        )

        # Create Standard Deviation Measurements for the primary phase EDS composition
        eds_measured_composition_stddev = Measurement(
            name="EDS Measured Composition StdDev",
            template=eds_stddev_measurement_template,
            parameters=[
                Parameter('Al', NominalReal(data['Results']['Measured Composition StdDev (%)']['Al'], '%'), template=eds_stddev_measurement_template),
                Parameter('Co', NominalReal(data['Results']['Measured Composition StdDev (%)']['Co'], '%'), template=eds_stddev_measurement_template),
                Parameter('Cr', NominalReal(data['Results']['Measured Composition StdDev (%)']['Cr'], '%'), template=eds_stddev_measurement_template),
                Parameter('Cu', NominalReal(data['Results']['Measured Composition StdDev (%)']['Cu'], '%'), template=eds_stddev_measurement_template),
                Parameter('Fe', NominalReal(data['Results']['Measured Composition StdDev (%)']['Fe'], '%'), template=eds_stddev_measurement_template),
                Parameter('Mn', NominalReal(data['Results']['Measured Composition StdDev (%)']['Mn'], '%'), template=eds_stddev_measurement_template),
                Parameter('Ni', NominalReal(data['Results']['Measured Composition StdDev (%)']['Ni'], '%'), template=eds_stddev_measurement_template),
                Parameter('V', NominalReal(data['Results']['Measured Composition StdDev (%)']['V'], '%'), template=eds_stddev_measurement_template)
            ]
        )

        # Create Measurements for 2nd Phase EDS composition results
        eds_2nd_phase_measured_composition = Measurement(
            name="2nd Phase EDS Measured Composition",
            template=eds_composition_measurement_template,  # Reuse the same template
            parameters=[
                Parameter('Al', NominalReal(data['Results']['2nd Phase EDS Measured Composition (%)']['Al'], '%'), template=eds_composition_measurement_template),
                Parameter('Co', NominalReal(data['Results']['2nd Phase EDS Measured Composition (%)']['Co'], '%'), template=eds_composition_measurement_template),
                Parameter('Cr', NominalReal(data['Results']['2nd Phase EDS Measured Composition (%)']['Cr'], '%'), template=eds_composition_measurement_template),
                Parameter('Cu', NominalReal(data['Results']['2nd Phase EDS Measured Composition (%)']['Cu'], '%'), template=eds_composition_measurement_template),
                Parameter('Fe', NominalReal(data['Results']['2nd Phase EDS Measured Composition (%)']['Fe'], '%'), template=eds_composition_measurement_template),
                Parameter('Mn', NominalReal(data['Results']['2nd Phase EDS Measured Composition (%)']['Mn'], '%'), template=eds_composition_measurement_template),
                Parameter('Ni', NominalReal(data['Results']['2nd Phase EDS Measured Composition (%)']['Ni'], '%'), template=eds_composition_measurement_template),
                Parameter('V', NominalReal(data['Results']['2nd Phase EDS Measured Composition (%)']['V'], '%'), template=eds_composition_measurement_template)
            ]
        )

        # Create Standard Deviation Measurements for the 2nd Phase EDS composition
        eds_2nd_phase_measured_composition_stddev = Measurement(
            name="2nd Phase EDS Measured Composition StdDev",
            template=eds_stddev_measurement_template,  # Reuse the same template for standard deviation
            parameters=[
                Parameter('Al', NominalReal(data['Results']['2nd Phase EDS Measured Composition StdDev (%)']['Al'], '%'), template=eds_stddev_measurement_template),
                Parameter('Co', NominalReal(data['Results']['2nd Phase EDS Measured Composition StdDev (%)']['Co'], '%'), template=eds_stddev_measurement_template),
                Parameter('Cr', NominalReal(data['Results']['2nd Phase EDS Measured Composition StdDev (%)']['Cr'], '%'), template=eds_stddev_measurement_template),
                Parameter('Cu', NominalReal(data['Results']['2nd Phase EDS Measured Composition StdDev (%)']['Cu'], '%'), template=eds_stddev_measurement_template),
                Parameter('Fe', NominalReal(data['Results']['2nd Phase EDS Measured Composition StdDev (%)']['Fe'], '%'), template=eds_stddev_measurement_template),
                Parameter('Mn', NominalReal(data['Results']['2nd Phase EDS Measured Composition StdDev (%)']['Mn'], '%'), template=eds_stddev_measurement_template),
                Parameter('Ni', NominalReal(data['Results']['2nd Phase EDS Measured Composition StdDev (%)']['Ni'], '%'), template=eds_stddev_measurement_template),
                Parameter('V', NominalReal(data['Results']['2nd Phase EDS Measured Composition StdDev (%)']['V'], '%'), template=eds_stddev_measurement_template)
            ]
        )

        # Create the overall sequence for EBSD and EDS Mapping Experiment
        ebsd_eds_mapping_sequence = MaterialsSequence(
            name="EBSD and EDS Mapping Experiment",
            science_kit=science_kit,
            material=ebsd_sample_material,
            ingredients=[sample_ingredient],
            process=ebsd_eds_mapping_process,
            measurements=[
                eds_measured_composition,
                eds_measured_composition_stddev,
                eds_2nd_phase_measured_composition,
                eds_2nd_phase_measured_composition_stddev
            ]
        )

        # Link internal elements within the sequence
        ebsd_eds_mapping_sequence.link_within()

        return ebsd_eds_mapping_sequence
    
    ebsd_eds_mapping_sequence = make_ebsd_eds_mapping_sequence(ebsd_eds_data)
    return ebsd_eds_mapping_sequence

class BIRDSHOTModeller(GEMDModeller):

    def __init__(self, files_folder, gemd_folder, instantiate_build):
        """
        Initialize the GEMDModeller with stores_config, files_folder, and gemd_folder.
        """
        super().__init__(files_folder, gemd_folder, instantiate_build)
        self.add_automatable_component(
            lambda s: 'NI-HSR' in s and (not ('.' in s)),  
            (r'\b[A-Z]{3}[0-9]{2}\b', True),
            [], 
            lambda file_name, file_path, component: ni_model(file_name, file_path, component)
        )
        self.add_automatable_component(
            lambda s: 'EDS' in s and (not ('.' in s)),  
            (r'\b[A-Z]{3}[0-9]{2}\b', True),
            [], 
            lambda file_name, file_path, component: eds_model(file_name, file_path, component)
        )
        self.start_folder_monitoring()

def main(args=None):
    """
    Main method to run from command line
    """
    BIRDSHOTModeller.run_from_command_line(args)

if __name__ == "__main__":
    main()
