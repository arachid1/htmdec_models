from openmsimodel.interactive.gemd_modeller import GEMDModeller
from gemd import (
    MaterialTemplate,
    ProcessTemplate,
    MeasurementTemplate,
    ParameterTemplate,
    ConditionTemplate,
    MaterialRun,
    MaterialSpec,
    RealBounds,
    CategoricalBounds,
)
from openmsimodel.structures.materials_sequence import MaterialsSequence
from openmsimodel.science_kit.science_kit import ScienceKit
from openmsimodel.entity.gemd.material import Material
from openmsimodel.entity.gemd.process import Process
from openmsimodel.entity.gemd.measurement import Measurement
from openmsimodel.entity.gemd.ingredient import Ingredient
import json
import os


def parse_control_script(file_name, file_path, component):
    """Parse a control script JSON file and construct a GEMD model for the machine operation."""

    science_kit = ScienceKit()

    # Load JSON data
    with open(file_path, "r") as f:
        data = json.load(f)

    # Extract relevant fields
    instructions = data.get("instructions", [])
    if not instructions:
        print(f"Skipping {file_name}: No instructions found.")
        return

    operation_type = instructions[0][0]  # e.g., "maxima"
    operation_settings = instructions[0][1]

    sample_data = operation_settings.get("sample", {})
    machine_settings = (
        operation_settings.get("xray_settings", {})
        or operation_settings.get("engraver_settings", {})
        or operation_settings.get("laser_settings", {})
    )
    collection_settings = operation_settings.get("collection_settings", {})

    # Define templates
    sample_material_template = MaterialTemplate("Sample")
    machine_process_template = ProcessTemplate(f"{operation_type} Process")
    machine_measurement_template = MeasurementTemplate(
        f"{operation_type} Measurement",
        parameters=[
            ParameterTemplate("Machine Current", bounds=RealBounds(0, 100, "mA")),
            ParameterTemplate("Machine Voltage", bounds=RealBounds(0, 500, "V")),
            ParameterTemplate("Beam X", bounds=RealBounds(0, 500, "mm")),
            ParameterTemplate("Beam Y", bounds=RealBounds(0, 500, "mm")),
        ],
    )

    # Create GEMD Objects
    sample_material = Material("Sample", template=sample_material_template)
    machine_process = Process(
        f"{operation_type} Execution", template=machine_process_template
    )
    machine_measurement = Measurement(
        f"{operation_type} Output", template=machine_measurement_template
    )

    # Assign parameter values
    machine_measurement.parameters = {
        "Machine Current": machine_settings.get("current", 0),
        "Machine Voltage": machine_settings.get("voltage", 0),
        "Beam X": machine_settings.get("beam_x", 0),
        "Beam Y": machine_settings.get("beam_y", 0),
    }

    # Create a sequence linking material, process, and measurement
    operation_sequence = MaterialsSequence(
        name=f"{operation_type} Control Sequence",
        science_kit=science_kit,
        ingredients=[Ingredient("Sample")],
        process=machine_process,
        material=sample_material,
        measurements=[machine_measurement],
    )
    operation_sequence.link_within()

    return science_kit.assets()


class ControlScriptModeller(GEMDModeller):

    def __init__(self, files_folder, gemd_folder, instantiate_build):
        """
        Initialize the GEMDModeller with stores_config, files_folder, and gemd_folder.
        """
        super().__init__(files_folder, gemd_folder, instantiate_build)

        # Add JSON control script parsing
        self.add_automatable_component(
            lambda s: s.endswith(".json"),
            (r".+\.json$", False),
            [],
            lambda file_name, file_path, component: parse_control_script(
                file_name, file_path, component
            ),
        )

        self.start_folder_monitoring()


def main(args=None):
    """
    Main method to run from the command line.
    """
    ControlScriptModeller.run_from_command_line(args)


if __name__ == "__main__":
    main()
