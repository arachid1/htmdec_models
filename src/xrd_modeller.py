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

def xrd_model(file_name, file_path, component):

    science_kit = ScienceKit()

    preprocessing_process_template = ProcessTemplate("Preprocessing")
    sample_material_template = MaterialTemplate("Sample")
    

    preprocessing_process = Process("Peprocessing", template=preprocessing_process_template)
    sample_material = Material("Sample", template=sample_material_template)
    preprocessed_sample_ingredient = Ingredient("Unprocessed Sample")

    external_preprocessing_sequence = MaterialsSequence(
        name=f"External Preprocessing",
        science_kit=science_kit,
        ingredients=[preprocessed_sample_ingredient],
        process=preprocessing_process,
        material=sample_material,
        measurements=[],
    )
    external_preprocessing_sequence.link_within()

    #####

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

    sample_ingredient_name = "Sample"
    sample_ingredient = Ingredient(sample_ingredient_name)
    thining_process = Process("Thining", template=thinning_process_template)
    thin_sample_material = Material("Thin Sample", template=sample_material_template)

    thining_sequence = MaterialsSequence(
        name=f"Thining Process",
        science_kit=science_kit,
        ingredients=[sample_ingredient],
        process=thining_process,
        material=thin_sample_material,
        measurements=[],
    )
    thining_sequence.link_within()
    thining_sequence.link_prior(external_preprocessing_sequence, ingredient_name_to_link=sample_ingredient_name)

    # ####

    polishing_process_template = ProcessTemplate("Polishing",
            parameters=ParameterTemplate(
                name="Grit Size",
                bounds=CategoricalBounds(['400', '600', '800'])
            ),
    )

    vibrometry_measurement_template = MeasurementTemplate('Vibrometry')
    edm_measurement_template = MeasurementTemplate('Electrical Discharge Machining')
    thin_sample_ingredient_name = "Thin Sample"
    thin_sample_ingredient = Ingredient(thin_sample_ingredient_name)
    polishing_process = Process("Polishing", template=polishing_process_template)
    polished_sample_material = Material("Polished Thin Sample", template=sample_material_template)
    vibrometry_measurement = Measurement("Vibrometry", template=vibrometry_measurement_template)
    edm_measurement = Measurement("Electrical Discharge Machining", template=edm_measurement_template)

    xrd_measurement_template = MeasurementTemplate("XRD", 
        parameters=[
            ParameterTemplate('Coordinate',
                description="X or Y coordinate of the sample based on a grid",
                bounds=RealBounds(0, 20, "")
            ),
            ParameterTemplate('Step Size',
                description="Step between two consecuctive exposures",
                bounds=RealBounds(0, 4, "")
            )]
    )
    xrd_measurement = Measurement("Single Shot Test", template=xrd_measurement_template)

    polishing_sequence = MaterialsSequence(
        name=f"Polishing Process",
        science_kit=science_kit,
        ingredients=[thin_sample_ingredient],
        process=polishing_process,
        material=polished_sample_material,
        measurements=[vibrometry_measurement, edm_measurement, xrd_measurement],
    )
    polishing_sequence.link_within()
    polishing_sequence.link_prior(thining_sequence, ingredient_name_to_link=thin_sample_ingredient_name)


    return science_kit.assets()

class XRDModeller(GEMDModeller):


    def __init__(self, files_folder, gemd_folder):
        """
        Initialize the GEMDModeller with stores_config, files_folder, and gemd_folder.
        """
        super().__init__(files_folder, gemd_folder)
        self.add_automatable_component(
            lambda s: True,  
            r'^([A-Za-z]+_\d+)',
            [], 
            lambda file_name, file_path, component:xrd_model(file_name, file_path, component)
        )


def main(args=None):
    """
    Main method to run from command line
    """
    XRDModeller.run_from_command_line(args)


if __name__ == "__main__":
    main()
