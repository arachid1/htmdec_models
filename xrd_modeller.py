from openmsimodel.interactive.gemd_modeller import GEMDModeller
from gemd import MaterialTemplate, ProcessTemplate, MeasurementTemplate, ParameterTemplate, MaterialRun, MaterialSpec, RealBounds
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

def xrd_model(file_name, component):
    preprocessing_process_template = ProcessTemplate("Preprocessing")
    sample_template = MaterialTemplate("Sample")
    science_kit = ScienceKit()

    preprocessing_process = Process("Peprocessing", template=preprocessing_process_template)
    sample_material = Material("Sample", template=sample_template)

    external_preprocessing_sequence = MaterialsSequence(
        name=f"External Preprocessing",
        science_kit=science_kit,
        ingredients=[],
        process=preprocessing_process,
        material=sample_material,
        measurements=[],
    )
    external_preprocessing_sequence.link_within()

    thinning_process_template = ProcessTemplate("Thinning", 
            parameters=ParameterTemplate(
                name="Tool",
                bounds=CategoricalBounds(["Diamond Wire Saw"])
            ),
            conditions=ConditionTemplate(
                name="Target Thickness",
                description="Objective of the thinning process which serves to achieve a target absorption length of the material",
                bounds=RealBounds(0, 20, "cm") #TODO 
        )
    )

    preprocessing_process = Process("Thining", template=thinning_process_template)


class XRDModeller(GEMDModeller):


    def __init__(self, **params):
        """
        Initialize the GEMDModeller with stores_config, files_folder, and gemd_folder.
        """
        super()._init__(**params)
        self.add_automatable_component(
            lambda s: s.endswith('.txt'),  # Rule function
            r"\d+",  # Pattern for extracting ID (you may want to modify this in the action function)
            [MaterialTemplate("Alloy")],  # Schema (or component)
            lambda file_name, component: Material(
                f'{self.extract_id_from_filename(component["pattern"], file_name)} Alloy', 
                template=component['schema'][0]  
            )
        )
        self.add_automatable_component(
            lambda s: s.endswith('.json'), 
            r"\d+",  
            [MeasurementTemplate("XRD")],  
            lambda file_name, component: Measurement(
                f'{self.extract_id_from_filename(component["pattern"], file_name)} XRD Analysis', 
                template=component['schema'][0]  
            )
        )


def main(args=None):
    """
    Main method to run from command line
    """
    XRDModeller.run_from_command_line(args)


if __name__ == "__main__":
    main()
