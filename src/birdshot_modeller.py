from openmsimodel.interactive.gemd_modeller import GEMDModeller
from openmsimodel.interactive.gemd_modeller import AutomatableComponentTree

from gemd import (
    MaterialTemplate,
    ProcessTemplate,
    MeasurementTemplate,
    ParameterTemplate,
    PropertyTemplate,
    ConditionTemplate,
    MaterialRun,
    MaterialSpec,
    RealBounds,
    CategoricalBounds,
    Parameter,
    NominalReal,
    FileLink,
    NominalCategorical,
    Condition,
    Property,
    PropertyAndConditions,
)
from openmsimodel.structures.materials_sequence import MaterialsSequence
from openmsimodel.science_kit.science_kit import ScienceKit
from openmsimodel.entity.gemd.material import Material
from openmsimodel.entity.gemd.process import Process
from openmsimodel.entity.gemd.measurement import Measurement
from openmsimodel.entity.gemd.ingredient import Ingredient
from openmsimodel.structures.materials_sequence import MaterialsSequence
import re
import os


class BIRDSHOTAutomatableComponentTree(AutomatableComponentTree):

    def __init__(self):
        super().__init__()
        self.schema = "dag"

    def map_to_schema(self):
        if self.schema == "dag":
            return self.map_to_dag()

    def map_to_dag(self):
        _all = []
        if "Syn" in self.file_mappings:
            synthesis_kit = self.file_mappings["Syn"].output
            if "EDS" in self.file_mappings:
                eds_kit = self.file_mappings["EDS"].output
                eds_kit_first_sequence = list(eds_kit.structures.values())[0]
                eds_kit.link_prior(
                    synthesis_kit,
                    ingredient_name_to_link=eds_kit_first_sequence.element_assets[
                        0
                    ].name,
                )
                _all.extend(eds_kit.assets())
                del eds_kit
            if "NI-HSR" in self.file_mappings:
                ni_hsr_kit = self.file_mappings["NI-HSR"].output
                ni_hsr_kit_first_sequence = list(ni_hsr_kit.structures.values())[0]
                ni_hsr_kit.link_prior(
                    synthesis_kit,
                    ingredient_name_to_link=ni_hsr_kit_first_sequence.element_assets[
                        0
                    ].name,
                )
                _all.extend(ni_hsr_kit.assets())
                del ni_hsr_kit
            _all.extend(synthesis_kit.assets())
        return _all


class BIRDSHOTModeller(GEMDModeller):

    def __init__(
        self,
        mode,
        files_folder,
        gemd_folder,
        api_url,
        girder_api_key,
        girder_root_folder_id,
        instantiate_build,
    ):
        """
        Initialize the GEMDModeller with stores_config, files_folder, and gemd_folder.
        """
        super().__init__(
            mode,
            files_folder,
            gemd_folder,
            api_url,
            girder_api_key,
            girder_root_folder_id,
            instantiate_build,
        )
        self.add_automatable_component(
            "Tensile",
            lambda file_name, file_path: "Tensile" in file_path,
            (r"\b[A-Z]{3}\d{2}_(?:VAM|DED)-[A-Z](?:_[A-Za-z]+_[a-z])?\b", False),
            lambda file_name, file_path, component, form: self.tensile_model(
                file_name, file_path, component, form
            ),
        )
        # self.add_automatable_component(
        #     "NI-HSR",
        #     lambda file_name, file_path: "NI-HSR" in file_path
        #     and (not ("." in file_path)),
        #     (r"\b[A-Z]{3}\d{2}_(?:VAM|DED)-[A-Z](?:_[A-Za-z]+_[a-z])?\b", False),
        #     lambda file_name, file_path, component, form: self.hsr_ni_model(
        #         file_name, file_path, component, form
        #     ),
        # )
        # self.add_automatable_component(
        #     "EDS",
        #     lambda file_name, file_path: "EDS" in file_path
        #     and (not ("." in file_path)),
        #     (r"\b[A-Z]{3}\d{2}_(?:VAM|DED)-[A-Z](?:_[A-Za-z]+_[a-z])?\b", False),
        #     lambda file_name, file_path, component, form: self.eds_model(
        #         file_name, file_path, component, form
        #     ),
        # )
        # self.add_automatable_component(
        #     "Syn",
        #     lambda file_name, file_path: "Syn" in file_path
        #     and (not ("." in file_path)),
        #     (r"\b[A-Z]{3}\d{2}_(?:VAM|DED)-[A-Z](?:_[A-Za-z]+_[a-z])?\b", False),
        #     lambda file_name, file_path, component, form: self.synthesis_model(
        #         file_name, file_path, component, form
        #     ),
        # )
        self.tree_class = BIRDSHOTAutomatableComponentTree
        self.start_monitoring()

    def tensile_model(self, file_name, file_path, component, form):
        science_kit = ScienceKit()
        tensile_data = form["data"]
        file_id = form["data"]["sampleId"]
        print(tensile_data)
        exit()

    def hsr_ni_model(self, file_name, file_path, component, form):
        science_kit = ScienceKit()
        ni_data = form["data"]
        file_id = form["data"]["sampleId"]

        files = [("File", item) for item in form["files"]]

        def make_sample_preparation_sequence():
            """Creates a Sample Preparation sequence (placeholder for now)."""

            sample_ingredient = Ingredient(f"{file_id} Sample for HSR-NI")

            sample_prep_process_template = ProcessTemplate("HSR-NI Sample Preparation")
            sample_prep_process = Process(
                f"{file_id} Sample Preparation Process",
                template=sample_prep_process_template,
            )

            prepared_sample_template = MaterialTemplate("HSR-NI Sample")
            prepared_sample = Material(
                f"{file_id} Prepared Sample", template=prepared_sample_template
            )

            # Create an empty sample preparation sequence
            sample_preparation_sequence = MaterialsSequence(
                name="Sample Preparation Sequence",
                science_kit=science_kit,
                ingredients=[sample_ingredient],
                process=sample_prep_process,
                material=prepared_sample,
                measurements=[],
            )

            # Link internal structure
            sample_preparation_sequence.link_within()

            return sample_preparation_sequence

        sample_preparation_sequence = make_sample_preparation_sequence()

        def make_nanoindentation_sequence(name):
            """Creates the Nanoindentation test sequence using metadata from Test & Analysis Parameters."""

            ni_sample_ingredient = Ingredient(name)

            # Extract metadata
            test_parameters = ni_data["Test Parameters"]
            analysis_parameters = ni_data["Analysis Parameters"]

            step_magnitute_parameter_template = ParameterTemplate(
                "Step Magnitude", bounds=RealBounds(0, 100, "mN")
            )
            target_load_parameter_template = ParameterTemplate(
                "Target Load", bounds=RealBounds(0, 20000, "mN")
            )
            ind_strain_rate_parameter_template = ParameterTemplate(
                "Target Ind Strain Rate", bounds=RealBounds(0, 10, "1/s")
            )
            load_time_parameter_template = ParameterTemplate(
                "Hold Maximum Load Time", bounds=RealBounds(0, 10, "seconds")
            )
            unload_time_parameter_template = ParameterTemplate(
                "Unload Time", bounds=RealBounds(0, 10, "seconds")
            )
            surface_stifness_trigger_template = ParameterTemplate(
                "Surface Stiffness Trigger", bounds=RealBounds(0, 1000, "N/m")
            )
            indent_number_parameter_template = ParameterTemplate(
                "Number of indents", bounds=RealBounds(0, 10, "dimensionless")
            )
            drift_parameter_template = ParameterTemplate(
                "Drift", bounds=RealBounds(0, 10, "nm/s")
            )
            drift_timeout_parameter_template = ParameterTemplate(
                "Drift Timeout", bounds=RealBounds(0, 1000, "seconds")
            )
            ####
            frame_stiffness_parameter_template = ParameterTemplate(
                "Frame Stiffness", bounds=RealBounds(0, 10**8, "N/m")
            )
            actuator_spring_stiffness_template = ParameterTemplate(
                "Actuator Spring Stiffness", bounds=RealBounds(0, 10**4, "N/m")
            )
            actuator_damping_coefficient_template = ParameterTemplate(
                "Actuator Damping Coefficient",
                bounds=RealBounds(0, 500, "N m/s"),
            )
            actuator_column_mass_template = ParameterTemplate(
                "Actuator Column Mass", bounds=RealBounds(0, 1, "kg")
            )
            actuator_calibration_coefficient_template = ParameterTemplate(
                "Actuator Calibration Coefficient",
                bounds=RealBounds(0, 100, "mV/mN"),
            )
            piezo_load_cell_stiffness_template = ParameterTemplate(
                "Piezo Load Cell Stiffness", bounds=RealBounds(0, 10**8, "N/m")
            )
            piezo_load_cell_damping_coefficient_template = ParameterTemplate(
                "Piezo Load Cell Damping Coefficient",
                bounds=RealBounds(0, 500, "N m/s"),
            )
            piezo_load_cell_effective_mass_template = ParameterTemplate(
                "Piezo Load Cell Effective Mass", bounds=RealBounds(0, 1, "kg")
            )
            piezo_calibration_coefficient_template = ParameterTemplate(
                "Piezo Calibration Coefficient",
                bounds=RealBounds(0, 10, "mV/mN"),
            )
            tip_area_coefficient_template = ParameterTemplate(
                "Tip Area Function Coefficient",
                bounds=RealBounds(0, 1, "dimensionless"),
            )

            # Initialize process template
            indentation_process_template = ProcessTemplate(
                "Nanoindentation Test",
                parameters=[
                    step_magnitute_parameter_template,
                    target_load_parameter_template,
                    ind_strain_rate_parameter_template,
                    load_time_parameter_template,
                    unload_time_parameter_template,
                    surface_stifness_trigger_template,
                    indent_number_parameter_template,
                    drift_parameter_template,
                    drift_timeout_parameter_template,
                    frame_stiffness_parameter_template,
                    actuator_spring_stiffness_template,
                    actuator_damping_coefficient_template,
                    actuator_column_mass_template,
                    actuator_calibration_coefficient_template,
                    piezo_load_cell_stiffness_template,
                    piezo_load_cell_damping_coefficient_template,
                    piezo_load_cell_effective_mass_template,
                    piezo_calibration_coefficient_template,
                    tip_area_coefficient_template,
                ],
            )
            indentation_process = Process(
                f"{file_id} Nanoindentation", template=indentation_process_template
            )

            indentation_process.update_parameters(
                Parameter(
                    "Step Magnitude",
                    value=NominalReal(test_parameters["Step Magnitude"], "mN"),
                    template=step_magnitute_parameter_template,
                ),
                which="run",
            )

            indentation_process.update_parameters(
                Parameter(
                    "Target Load",
                    value=NominalReal(test_parameters["Target Load"], "mN"),
                    template=target_load_parameter_template,
                ),
                which="run",
            )

            indentation_process.update_parameters(
                Parameter(
                    "Target Ind Strain Rate",
                    value=NominalReal(test_parameters["Target Ind Strain Rate"], "1/s"),
                    template=ind_strain_rate_parameter_template,
                ),
                which="run",
            )

            indentation_process.update_parameters(
                Parameter(
                    "Hold Maximum Load Time",
                    value=NominalReal(
                        test_parameters["Hold Maximum Load Time"], "seconds"
                    ),
                    template=load_time_parameter_template,
                ),
                which="run",
            )
            indentation_process.update_parameters(
                Parameter(
                    "Unload Time",
                    value=NominalReal(test_parameters["Unload Time"], "seconds"),
                    template=unload_time_parameter_template,
                ),
                which="run",
            )
            indentation_process.update_parameters(
                Parameter(
                    "Surface Stiffness Trigger",
                    value=NominalReal(
                        test_parameters["Surface Stiffness Trigger"], "N/s"
                    ),
                    template=surface_stifness_trigger_template,
                ),
                which="run",
            )
            indentation_process.update_parameters(
                Parameter(
                    "Number of indents",
                    value=NominalReal(
                        test_parameters["Number of Indents"], "dimensionless"
                    ),
                    template=indent_number_parameter_template,
                ),
                which="run",
            )

            indentation_process.update_parameters(
                Parameter(
                    "Drift",
                    value=NominalReal(test_parameters["Engage Options, Drift"], "nm/s"),
                    template=drift_parameter_template,
                ),
                which="run",
            )

            indentation_process.update_parameters(
                Parameter(
                    "Drift Timeout",
                    value=NominalReal(
                        test_parameters["Hold Maximum Load Time"], "seconds"
                    ),
                    template=load_time_parameter_template,
                ),
                which="run",
            )

            ###

            indentation_process.update_parameters(
                Parameter(
                    "Frame Stiffness",
                    value=NominalReal(analysis_parameters["Frame Stiffness"], "N/m"),
                    template=frame_stiffness_parameter_template,
                ),
                which="run",
            )

            indentation_process.update_parameters(
                Parameter(
                    "Actuator Spring Stiffness",
                    value=NominalReal(
                        analysis_parameters["Actuator Spring Stiffness"], "N/m"
                    ),
                    template=actuator_spring_stiffness_template,
                ),
                which="run",
            )

            indentation_process.update_parameters(
                Parameter(
                    "Actuator Damping Coefficient",
                    value=NominalReal(
                        analysis_parameters["Actuator Damping Coefficient"], "N m/s"
                    ),
                    template=actuator_damping_coefficient_template,
                ),
                which="run",
            )

            indentation_process.update_parameters(
                Parameter(
                    "Actuator Column Mass",
                    value=NominalReal(
                        analysis_parameters["Actuator Column Mass"], "kg"
                    ),
                    template=actuator_column_mass_template,
                ),
                which="run",
            )

            indentation_process.update_parameters(
                Parameter(
                    "Actuator Calibration Coefficient",
                    value=NominalReal(
                        analysis_parameters["Actuator Calibration Coefficient"],
                        "dimensionless",
                    ),
                    template=actuator_calibration_coefficient_template,
                ),
                which="run",
            )

            indentation_process.update_parameters(
                Parameter(
                    "Piezo Load Cell Stiffness",
                    value=NominalReal(
                        analysis_parameters["Pieze Load Cell Stiffness"], "N/m"
                    ),
                    template=piezo_load_cell_stiffness_template,
                ),
                which="run",
            )

            indentation_process.update_parameters(
                Parameter(
                    "Piezo Load Cell Damping Coefficient",
                    value=NominalReal(
                        analysis_parameters["Piezo Load Cell Damping Coefficient"],
                        "N m/s",
                    ),
                    template=piezo_load_cell_damping_coefficient_template,
                ),
                which="run",
            )

            indentation_process.update_parameters(
                Parameter(
                    "Piezo Load Cell Effective Mass",
                    value=NominalReal(
                        analysis_parameters["Piezo Load Cell Effective Mass"], "kg"
                    ),
                    template=piezo_load_cell_effective_mass_template,
                ),
                which="run",
            )

            indentation_process.update_parameters(
                Parameter(
                    "Piezo Calibration Coefficient",
                    value=NominalReal(
                        analysis_parameters["Piezo Calibration Coefficient"],
                        "dimensionless",
                    ),
                    template=piezo_calibration_coefficient_template,
                ),
                which="run",
            )
            for coeff, value in analysis_parameters[
                "Tip Area Function Coefficients"
            ].items():
                indentation_process.update_parameters(
                    Parameter(
                        f"Tip Area Function {coeff}",
                        value=NominalReal(value, "dimensionless"),
                        template=tip_area_coefficient_template,  # Single shared template
                    ),
                    which="run",
                )

            nanoindentation_sequence = MaterialsSequence(
                name="Nanoindentation Test Sequence",
                science_kit=science_kit,
                ingredients=[ni_sample_ingredient],
                process=indentation_process,
                measurements=[],
            )

            # Link internal structure
            nanoindentation_sequence.link_within()

            return nanoindentation_sequence

        ing_name = f"{file_id} Prepared Sample for HSR-NI"
        nanoindentation_sequence = make_nanoindentation_sequence(ing_name)
        nanoindentation_sequence.link_prior(
            sample_preparation_sequence, ingredient_name_to_link=ing_name
        )
        return science_kit

    def synthesis_model(self, file_name, file_path, component, form):

        science_kit = ScienceKit()
        if not form:
            match = re.search(component["file_id_regex_pattern"][0], file_path)

            file_id = match.group()
            if not match:
                print("No pattern found.")
                return

            raw_data_forms = self.client.get(
                "entry/search",
                parameters={"query": f"^{match.group()[:3]}.._VAM-.", "limit": 1000},
            )
            synthesis_data = {}
            for form in raw_data_forms:
                if (
                    "Syn" in form["data"]["targetPath"]
                    and match.group() in form["data"]["targetPath"]
                ):
                    synthesis_data.update(form["data"])
                    file_id = form["data"]["sampleId"]

            synthesis_folder_id = self.find_folder_by_path(
                self.girder_root_folder_id, form["data"]["targetPath"]
            )

            dms_syn_folder_items = self.client.get(
                f"/item",
                parameters={
                    "folderId": synthesis_folder_id,
                },
            )

            files = [(item["name"], item["_id"]) for item in dms_syn_folder_items]
        else:
            synthesis_data = form["data"]
            file_id = form["data"]["sampleId"]

            files = [(item["name"], item["_id"]) for item in form["files"]]

        def make_forging_sequence(data):
            ingot_ingredient = Ingredient(f"{file_id} Ingot")

            forging_process_template = ProcessTemplate("Forging")
            forging_process = Process(
                f"{file_id} Forging", template=forging_process_template
            )

            soak_time_property_template = PropertyTemplate(
                "Soak Time", bounds=RealBounds(0, 30, "minute")
            )
            temperature_property_template = PropertyTemplate(
                "Temperature", bounds=RealBounds(0, 30, "celsius")
            )
            ingot_material_template = MaterialTemplate(
                "Ingot",
                properties=[soak_time_property_template, temperature_property_template],
            )

            forged_ingot_material = Material(
                f"{file_id} Forged Ingot", template=ingot_material_template
            )

            # TODO: change to a measurement

            dimension_property_template = PropertyTemplate(
                "Dimension", bounds=RealBounds(0, 15, "cm")
            )
            dimension_measurement_template = MeasurementTemplate(
                "Dimensions", properties=dimension_property_template
            )
            dimensions_before_measurement = Measurement(
                f"{file_id} Prior Dimension", template=dimension_measurement_template
            )

            if "Forging" in data:
                try:
                    dimensions_before_measurement.update_properties(  # TODO: add thickness reduction
                        Property(
                            "Prior Length",
                            value=NominalReal(
                                data["Forging"]["Ingot Dimensions Before"]["Length"],
                                "cm",
                            ),
                            template=dimension_property_template,
                        ),
                        which="run",
                    )
                except Exception as e:
                    print(f"Error processing attribute: {e}. Skipping...")
                try:
                    dimensions_before_measurement.update_properties(  # TODO: add thickness reduction
                        Property(
                            "Prior Thickness",
                            value=NominalReal(
                                synthesis_data["Forging"]["Ingot Dimensions Before"][
                                    "Thickness"
                                ],
                                "cm",
                            ),
                            template=dimension_property_template,
                        ),
                        which="run",
                    )
                except Exception as e:
                    print(f"Error processing attribute: {e}. Skipping...")
                try:
                    dimensions_before_measurement.update_properties(  # TODO: add thickness reduction
                        Property(
                            "Prior Width",
                            value=NominalReal(
                                synthesis_data["Forging"]["Ingot Dimensions Before"][
                                    "Width"
                                ],
                                "cm",
                            ),
                            template=dimension_property_template,
                        ),
                        which="run",
                    )
                except Exception as e:
                    print(f"Error processing attribute: {e}. Skipping...")
            dimensions_after_measurement = Measurement(
                f"{file_id} Posterior Dimension",
                template=dimension_measurement_template,
            )
            if "Forging" in synthesis_data:
                try:
                    dimensions_after_measurement.update_properties(  # TODO: add thickness reduction
                        Property(
                            "Posterior Length",
                            value=NominalReal(
                                synthesis_data["Forging"]["Ingot Dimensions After"][
                                    "Length"
                                ],
                                "cm",
                            ),
                            template=dimension_property_template,
                        ),
                    )
                except Exception as e:
                    print(f"Error processing attribute: {e}. Skipping...")
                try:
                    dimensions_after_measurement.update_properties(  # TODO: add thickness reduction
                        Property(
                            "Posterior Length",
                            value=NominalReal(
                                synthesis_data["Forging"]["Ingot Dimensions After"][
                                    "Length"
                                ],
                                "cm",
                            ),
                            template=dimension_property_template,
                        ),
                    )
                except Exception as e:
                    print(f"Error processing attribute: {e}. Skipping...")
                try:
                    dimensions_after_measurement.update_properties(  # TODO: add thickness reduction
                        Property(
                            "Posterior Thickness",
                            value=NominalReal(
                                synthesis_data["Forging"]["Ingot Dimensions After"][
                                    "Thickness"
                                ],
                                "cm",
                            ),
                            template=dimension_property_template,
                        ),
                    )
                except Exception as e:
                    print(f"Error processing attribute: {e}. Skipping...")
                try:
                    dimensions_after_measurement.update_properties(  # TODO: add thickness reduction
                        Property(
                            "Posterior Width",
                            value=NominalReal(
                                synthesis_data["Forging"]["Ingot Dimensions After"][
                                    "Width"
                                ],
                                "cm",
                            ),
                            template=dimension_property_template,
                        ),
                        which="run",
                    )
                except Exception as e:
                    print(f"Error processing attribute: {e}. Skipping...")

            soak_time_property_template = PropertyTemplate(
                "Soak Time", bounds=RealBounds(0, 30, "minute")
            )
            soak_time_measurement_template = MeasurementTemplate(
                "Soak Time", properties=soak_time_property_template
            )
            soak_time_measurement = Measurement(
                f"{file_id} Soak Time", template=soak_time_measurement_template
            )

            if "Forging" in data:
                try:
                    soak_time_measurement.update_properties(  # TODO: update_properties causes clash and shouldnt be called
                        Property(
                            "Soak Time",
                            value=NominalReal(
                                data["Forging"]["Ingot Condition"]["Soak Time"],
                                "minute",
                            ),
                            template=soak_time_property_template,
                        ),
                        which="run",
                    )
                except Exception as e:
                    print(f"Error processing attribute: {e}. Skipping...")

            temperature_property_template = PropertyTemplate(
                "Temperature", bounds=RealBounds(0, 30, "celsius")
            )
            temperature_measurement_template = MeasurementTemplate(
                "Temperature", properties=temperature_property_template
            )
            temperature_measurement = Measurement(
                f"{file_id} Temperature", template=temperature_measurement_template
            )

            if "Forging" in data:
                try:
                    temperature_measurement.update_properties(  # TODO: update_properties causes clash and shouldnt be called
                        Property(
                            "Temperature",
                            value=NominalReal(
                                synthesis_data["Forging"]["Ingot Condition"][
                                    "Temperature"
                                ],
                                "celsius",
                            ),
                            template=temperature_property_template,
                        ),
                        which="run",
                    )
                except Exception as e:

                    print(f"Error processing attribute: {e}. Skipping...")
                    if file_id == "AAA15_VAM-B":
                        print(synthesis_data["Forging"]["Ingot Condition"])
                        exit()

            forging_sequence = MaterialsSequence(
                name="Ingot Forging Sequence",
                science_kit=science_kit,
                ingredients=[ingot_ingredient],
                material=forged_ingot_material,
                process=forging_process,
                measurements=[
                    dimensions_before_measurement,
                    dimensions_after_measurement,
                    soak_time_measurement,
                    temperature_measurement,
                ],
            )

            forging_sequence.link_within()
            return forging_sequence

        forging_sequence = make_forging_sequence(synthesis_data)

        def make_homogenization_sequence(data, ingredient_name):

            forged_sample_ingredient = Ingredient(ingredient_name)
            atmosphere_template = ParameterTemplate(
                name="Atmosphere", bounds=CategoricalBounds(["Ar", "O2", "N2"])
            )

            cooling_rate_template = ParameterTemplate(
                name="Cooling Rate",
                bounds=CategoricalBounds(
                    ["FC", "AC", "WC"]
                ),  # FC: Furnace Cooling, AC: Air Cooling, WC: Water Cooling
            )

            duration_template = ParameterTemplate(
                name="Duration", bounds=RealBounds(0, 100, "hours")
            )

            pressure_template = ParameterTemplate(
                name="Pressure", bounds=RealBounds(0, 1000, "Pa")
            )

            temperature_template = ParameterTemplate(
                name="Temperature", bounds=RealBounds(0, 2000, "Celsius")
            )

            purging_pressure_template = ParameterTemplate(
                name="Purging Sequence Pressure", bounds=RealBounds(0, 100, "Pa")
            )

            # Create Homogenization process template
            homogenization_process_template = ProcessTemplate(
                name="Homogenization",
                parameters=[
                    atmosphere_template,
                    cooling_rate_template,
                    duration_template,
                    pressure_template,
                    temperature_template,
                    purging_pressure_template,
                ],
            )

            # Create Homogenization process
            homogenization_process = Process(
                name=f"{file_id} Homogenization",
                template=homogenization_process_template,
            )

            # Update the parameters with values from the provided data
            if "Homogenization" in data:
                try:
                    homogenization_process.update_parameters(
                        Parameter(
                            "Atmosphere",
                            value=NominalCategorical(
                                data["Homogenization"]["Thermal Conditions"][
                                    "Atmosphere"
                                ]
                            ),
                            template=atmosphere_template,
                        ),
                        which="run",
                    )
                except Exception as e:
                    print(f"Error processing attribute: {e}. Skipping...")
                try:
                    homogenization_process.update_parameters(
                        Parameter(
                            "Cooling Rate",
                            value=NominalCategorical(
                                data["Homogenization"]["Thermal Conditions"][
                                    "Cooling Rate"
                                ]
                            ),
                            template=cooling_rate_template,
                        ),
                        which="run",
                    )
                except Exception as e:
                    print(f"Error processing attribute: {e}. Skipping...")
                try:
                    homogenization_process.update_parameters(
                        Parameter(
                            "Duration",
                            value=NominalReal(
                                data["Homogenization"]["Thermal Conditions"][
                                    "Duration"
                                ],
                                "hours",
                            ),
                            template=duration_template,
                        ),
                        which="run",
                    )
                except Exception as e:
                    print(f"Error processing attribute: {e}. Skipping...")
                try:
                    homogenization_process.update_parameters(
                        Parameter(
                            "Pressure",
                            value=NominalReal(
                                data["Homogenization"]["Thermal Conditions"][
                                    "Pressure"
                                ],
                                "Pa",
                            ),
                            template=pressure_template,
                        ),
                        which="run",
                    )
                except Exception as e:
                    print(f"Error processing attribute: {e}. Skipping...")
                try:
                    homogenization_process.update_parameters(
                        Parameter(
                            "Pressure",
                            value=NominalReal(
                                data["Homogenization"]["Thermal Conditions"][
                                    "Pressure"
                                ],
                                "Pa",
                            ),
                            template=pressure_template,
                        ),
                        which="run",
                    )
                except Exception as e:
                    print(f"Error processing attribute: {e}. Skipping...")
                try:
                    homogenization_process.update_parameters(
                        Parameter(
                            "Temperature",
                            value=NominalReal(
                                data["Homogenization"]["Thermal Conditions"][
                                    "Temperature"
                                ],
                                "Celsius",
                            ),
                            template=temperature_template,
                        ),
                        which="run",
                    )
                except Exception as e:
                    print(f"Error processing attribute: {e}. Skipping...")

                # Add purging sequence pressures
                if "Purging Sequence Pressure" in data["Homogenization"]:
                    for step, pressure in data["Homogenization"][
                        "Purging Sequence Pressure"
                    ].items():
                        try:
                            homogenization_process.update_parameters(
                                Parameter(
                                    f"Purging Sequence Pressure {step}",
                                    value=NominalReal(pressure, "Pa"),
                                    template=purging_pressure_template,
                                )
                            )
                        except Exception as e:
                            print(f"Error processing attribute: {e}. Skipping...")

            homogenenous_material = Material(
                f"{file_id} Homogenous Material", template=MaterialTemplate("Sample")
            )

            time_parameter_template = ParameterTemplate(
                name="Time", bounds=RealBounds(0, 24, "hours")
            )
            # Measurement for time spent on the process
            time_spent_measurement_template = MeasurementTemplate(
                name="Time Spent", parameters=time_parameter_template
            )

            time_spent_measurement = Measurement(
                name=f"{file_id} Time Spent Measurement",
                template=time_spent_measurement_template,
            )

            try:
                time_spent_measurement.update_parameters(
                    Parameter(
                        "Time",
                        value=NominalReal(5, "hours"),
                        template=time_parameter_template,
                    )
                )
            except Exception as e:
                print(f"Error processing attribute: {e}. Skipping...")

            homogenization_sequence = MaterialsSequence(
                name="Homogenization Process",
                science_kit=science_kit,
                ingredients=[forged_sample_ingredient],
                material=homogenenous_material,  # Specify material if necessary
                process=homogenization_process,
                measurements=[time_spent_measurement],
            )

            homogenization_sequence.link_within()
            return homogenization_sequence

        ingredient_name = f"{file_id} Forged Sample"
        homogenization_sequence = make_homogenization_sequence(
            synthesis_data, ingredient_name
        )
        homogenization_sequence.link_prior(
            forging_sequence, ingredient_name_to_link=ingredient_name
        )

        def make_arc_melting_sequence(data, ingredient_name):

            homogenous_sample_ingredient = Ingredient(ingredient_name)

            argon_pressure_template = ConditionTemplate(
                name="Argon Pressure", bounds=RealBounds(0, 1000, "pascal")
            )
            vacuum_before_melt_template = ConditionTemplate(
                name="Vacuum Before Melt", bounds=RealBounds(0, 1, "pascal")
            )
            arc_melting_process_template = ProcessTemplate(
                name="Arc Melting",
                conditions=[argon_pressure_template, vacuum_before_melt_template],
            )
            arc_melting_process = Process(
                name=f"{file_id} Arc Melting", template=arc_melting_process_template
            )
            for file in files:
                arc_melting_process.update_filelinks(
                    FileLink(
                        f"{file[0]} (DMS)",
                        url=f"https://data.htmdec.org/api/v1/item/{file[1]}",
                    ),
                    which="run",
                )

            try:
                arc_melting_process.update_conditions(
                    Condition(
                        "Argon Pressure",
                        value=NominalReal(
                            data["Arc Melting"]["VAM Details"]["Argon Pressure"],
                            "pascal",
                        ),
                        template=argon_pressure_template,
                    ),
                    which="run",
                )
            except Exception as e:
                print(f"Error processing attribute: {e}. Skipping...")
            try:
                arc_melting_process.update_conditions(
                    Condition(
                        "Vacuum Before Melt",
                        value=NominalReal(
                            data["Arc Melting"]["VAM Details"]["Vacuum Before Melt"],
                            "pascal",
                        ),
                        template=vacuum_before_melt_template,
                    ),
                    which="run",
                )
            except Exception as e:
                print(f"Error processing attribute: {e}. Skipping...")

            mass_property_template = PropertyTemplate(
                name="Mass in grams",
                bounds=RealBounds(0, 100, "gram"),  # RealBounds for mass in grams
            )
            sample_material_template = MaterialTemplate(
                "Sample", properties=mass_property_template
            )
            arc_melted_material = Material(
                name=f"{file_id} Arc Melted Sample", template=sample_material_template
            )

            try:
                arc_melted_material.update_properties_and_conditions(  # TODO: update_properties causes clash and shouldnt be called
                    PropertyAndConditions(
                        property=Property(
                            "Target Al Mass",
                            value=NominalReal(
                                data["Material Preparation"]["Target Mass"]["Al"],
                                "gram",
                            ),
                            template=mass_property_template,
                        ),
                        conditions=[],
                    ),
                )
            except Exception as e:
                print(f"Error processing attribute: {e}. Skipping...")
            try:
                arc_melted_material.update_properties_and_conditions(  # TODO: update_properties causes clash and shouldnt be called
                    PropertyAndConditions(
                        property=Property(
                            "Target Co Mass",
                            value=NominalReal(
                                data["Material Preparation"]["Target Mass"]["Co"],
                                "gram",
                            ),
                            template=mass_property_template,
                        ),
                        conditions=[],
                    ),
                )
            except Exception as e:
                print(f"Error processing attribute: {e}. Skipping...")
            try:
                arc_melted_material.update_properties_and_conditions(  # TODO: update_properties causes clash and shouldnt be called
                    PropertyAndConditions(
                        property=Property(
                            "Target Cr Mass",
                            value=NominalReal(
                                data["Material Preparation"]["Target Mass"]["Cr"],
                                "gram",
                            ),
                            template=mass_property_template,
                        ),
                        conditions=[],
                    ),
                )
            except Exception as e:
                print(f"Error processing attribute: {e}. Skipping...")
            try:
                arc_melted_material.update_properties_and_conditions(  # TODO: update_properties causes clash and shouldnt be called
                    PropertyAndConditions(
                        property=Property(
                            "Target Fe Mass",
                            value=NominalReal(
                                data["Material Preparation"]["Target Mass"]["Fe"],
                                "gram",
                            ),
                            template=mass_property_template,
                        ),
                        conditions=[],
                    ),
                )
            except Exception as e:
                print(f"Error processing attribute: {e}. Skipping...")
            try:
                arc_melted_material.update_properties_and_conditions(  # TODO: update_properties causes clash and shouldnt be called
                    PropertyAndConditions(
                        property=Property(
                            "Target Mn Mass",
                            value=NominalReal(
                                data["Material Preparation"]["Target Mass"]["Mn"],
                                "gram",
                            ),
                            template=mass_property_template,
                        ),
                        conditions=[],
                    ),
                )
            except Exception as e:
                print(f"Error processing attribute: {e}. Skipping...")
            try:
                arc_melted_material.update_properties_and_conditions(  # TODO: update_properties causes clash and shouldnt be called
                    PropertyAndConditions(
                        property=Property(
                            "Target Ni Mass",
                            value=NominalReal(
                                data["Material Preparation"]["Target Mass"]["Ni"],
                                "gram",
                            ),
                            template=mass_property_template,
                        ),
                        conditions=[],
                    ),
                )
            except Exception as e:
                print(f"Error processing attribute: {e}. Skipping...")
            try:
                arc_melted_material.update_properties_and_conditions(  # TODO: update_properties causes clash and shouldnt be called
                    PropertyAndConditions(
                        property=Property(
                            "Target V Mass",
                            value=NominalReal(
                                data["Material Preparation"]["Target Mass"]["V"], "gram"
                            ),
                            template=mass_property_template,
                        ),
                        conditions=[],
                    ),
                )
            except Exception as e:
                print(f"Error processing attribute: {e}. Skipping...")

            # Create a measurement for the weighed mass
            weighed_mass_measurement_template = MeasurementTemplate(name="Weighing")
            weighed_mass_measurement = Measurement(
                name=f"{file_id} Weighed Mass",
                template=weighed_mass_measurement_template,
            )

            try:
                weighed_mass_measurement.update_properties(
                    Property(
                        "Al",
                        value=NominalReal(
                            data["Material Preparation"]["Weighed Mass"]["Al"], "gram"
                        ),
                        template=mass_property_template,
                    ),
                    which="run",
                )
            except Exception as e:
                print(f"Error processing attribute: {e}. Skipping...")
            try:
                weighed_mass_measurement.update_properties(
                    Property(
                        "Co",
                        value=NominalReal(
                            data["Material Preparation"]["Weighed Mass"]["Co"], "gram"
                        ),
                        template=mass_property_template,
                    ),
                    which="run",
                )
            except Exception as e:
                print(f"Error processing attribute: {e}. Skipping...")
            try:
                weighed_mass_measurement.update_properties(
                    Property(
                        "Cr",
                        value=NominalReal(
                            data["Material Preparation"]["Weighed Mass"]["Cr"], "gram"
                        ),
                        template=mass_property_template,
                    ),
                    which="run",
                )
            except Exception as e:
                print(f"Error processing attribute: {e}. Skipping...")
            try:
                weighed_mass_measurement.update_properties(
                    Property(
                        "Fe",
                        value=NominalReal(
                            data["Material Preparation"]["Weighed Mass"]["Fe"], "gram"
                        ),
                        template=mass_property_template,
                    ),
                    which="run",
                )
            except Exception as e:
                print(f"Error processing attribute: {e}. Skipping...")
            try:
                weighed_mass_measurement.update_properties(
                    Property(
                        "Mn",
                        value=NominalReal(
                            data["Material Preparation"]["Weighed Mass"]["Mn"], "gram"
                        ),
                        template=mass_property_template,
                    ),
                    which="run",
                )
            except Exception as e:
                print(f"Error processing attribute: {e}. Skipping...")
            try:
                weighed_mass_measurement.update_properties(
                    Property(
                        "Ni",
                        value=NominalReal(
                            data["Material Preparation"]["Weighed Mass"]["Ni"], "gram"
                        ),
                        template=mass_property_template,
                    ),
                    which="run",
                )
            except Exception as e:
                print(f"Error processing attribute: {e}. Skipping...")
            try:
                weighed_mass_measurement.update_properties(
                    Property(
                        "V",
                        value=NominalReal(
                            data["Material Preparation"]["Weighed Mass"]["V"], "gram"
                        ),
                        template=mass_property_template,
                    ),
                    which="run",
                )
            except Exception as e:
                print(f"Error processing attribute: {e}. Skipping...")

            # Create the MaterialsSequence for Arc Melting
            arc_melting_sequence = MaterialsSequence(
                name="Arc Melting Sequence",
                science_kit=science_kit,
                ingredients=[homogenous_sample_ingredient],
                material=arc_melted_material,
                process=arc_melting_process,
                measurements=[weighed_mass_measurement],
            )

            arc_melting_sequence.link_within()
            return arc_melting_sequence

        ingredient_name = f"{file_id} Homogenous Sample"
        arc_melting_sequence = make_arc_melting_sequence(
            synthesis_data, ingredient_name
        )
        arc_melting_sequence.link_prior(
            homogenization_sequence, ingredient_name_to_link=ingredient_name
        )

        return science_kit

    def eds_model(self, file_name, file_path, component, form=None):
        science_kit = ScienceKit()
        if not form:
            match = re.search(component["file_id_regex_pattern"][0], file_path)

            if not match:
                print("No pattern found.")
                return

            file_id = match.group()

            raw_data_forms = self.client.get(
                "entry/search",
                parameters={"query": f"^{file_id[:3]}.._VAM-.", "limit": 1000},
            )
            for form in raw_data_forms:
                if (
                    "EDS" in form["data"]["targetPath"]
                    and file_id in form["data"]["targetPath"]
                ):
                    ebsd_eds_data = form["data"]
                    file_id = form["data"]["sampleId"]

            eds_folder_id = self.find_folder_by_path(
                self.girder_root_folder_id, ebsd_eds_data["targetPath"]
            )

            dms_eds_folder_items = self.client.get(
                f"/item",
                parameters={
                    "folderId": eds_folder_id,
                },
            )

            files = [(item["name"], item["_id"]) for item in dms_eds_folder_items]
        else:
            ebsd_eds_data = form["data"]
            file_id = form["data"]["sampleId"]
            files = [(item["name"], item["_id"]) for item in form["files"]]

        def make_ebsd_eds_mapping_sequence(data):

            eds_measurements = []

            # Define material template for the sample
            sample_material_template = MaterialTemplate("Sample")

            # Ingredients (e.g., the material being analyzed)
            sample_ingredient = Ingredient(f"{file_id} EBS Ingredient")

            # Define individual Parameter Templates for the EBSD and EDS Mapping process
            beam_current_template = ParameterTemplate(
                name="Beam Current", bounds=RealBounds(-1, 100, "nanoampere")
            )
            beam_voltage_template = ParameterTemplate(
                name="Beam Voltage", bounds=RealBounds(0, 30, "kilovolt")
            )
            dwell_time_template = ParameterTemplate(
                name="Dwell Time", bounds=RealBounds(-1, 10, "seconds")
            )
            sample_tilt_template = ParameterTemplate(
                name="Sample Tilt", bounds=RealBounds(-1, 90, "degrees")
            )
            working_distance_template = ParameterTemplate(
                name="Working Distance", bounds=RealBounds(-1, 20, "millimeter")
            )
            low_vacuum_template = ParameterTemplate(
                name="Low Vacuum",
                bounds=CategoricalBounds(["None", "Low", "Medium", "High"]),
            )

            # Create process template for EBSD and EDS Mapping using the individual parameter templates
            ebsd_eds_mapping_process_template = ProcessTemplate(
                name="EBSD and EDS Mapping",
                parameters=[
                    beam_current_template,
                    beam_voltage_template,
                    dwell_time_template,
                    sample_tilt_template,
                    working_distance_template,
                    low_vacuum_template,
                ],
            )

            # Create EBSD and EDS Mapping process
            ebsd_eds_mapping_process = Process(
                name=f"{file_id} EBSD and EDS Mapping",
                template=ebsd_eds_mapping_process_template,
            )
            for file in files:
                ebsd_eds_mapping_process.update_filelinks(
                    FileLink(
                        f"{file[0]} (DMS)",
                        url=f"https://data.htmdec.org/api/v1/item/{file[1]}",
                    ),
                    which="run",
                )

            # Fill EBSD and EDS Mapping process with values from the form
            ebsd_eds_mapping_process.update_parameters(
                Parameter(
                    "Beam Current",
                    value=NominalReal(
                        data["EBSD and EDS Mapping"]["Beam Current"], "nanoampere"
                    ),
                    template=beam_current_template,
                ),
                Parameter(
                    "Beam Voltage",
                    value=NominalReal(
                        data["EBSD and EDS Mapping"]["Beam Voltage"], "kilovolt"
                    ),
                    template=beam_voltage_template,
                ),
                Parameter(
                    "Sample Tilt",
                    value=NominalReal(
                        data["EBSD and EDS Mapping"]["Sample Tilt"], "degrees"
                    ),
                    template=sample_tilt_template,
                ),
                Parameter(
                    "Working Distance",
                    value=NominalReal(
                        data["EBSD and EDS Mapping"]["Working Distance"], "millimeter"
                    ),
                    template=working_distance_template,
                ),
                Parameter(
                    "Low Vacuum",
                    value=NominalCategorical(
                        data["EBSD and EDS Mapping"]["Low Vacuum"]
                    ),
                    template=low_vacuum_template,
                ),
                which="run",
            )

            try:
                ebsd_eds_mapping_process.update_parameters(
                    Parameter(
                        "Dwell Time",
                        value=NominalReal(
                            data["EBSD and EDS Mapping"]["Dwell Time"], "seconds"
                        ),
                        template=dwell_time_template,
                    ),
                )
            except Exception as e:
                print(f"Error processing attribute: {e}. Skipping...")

            # Create Material for the EBSD and EDS Mapping sample
            ebsd_sample_material = Material(
                f"{file_id} Sample from EBSD and EDS Mapping",
                template=sample_material_template,
            )

            ### Create Measurement templates for EDS Measured Composition and StdDev
            eds_composition_measurement_template = MeasurementTemplate(
                "EDS Measured Composition"
            )
            eds_stddev_measurement_template = MeasurementTemplate(
                "EDS Measured Composition StdDev"
            )
            composition_parameter_template = ParameterTemplate(
                "Composition", bounds=RealBounds(0, 100, "")
            )
            composition_std_parameter_template = ParameterTemplate(
                "Composition Standard Deviation", bounds=RealBounds(0, 1, "")
            )

            # # Create Measurements for EDS composition results (primary phase)
            eds_measured_composition = Measurement(
                name=f"{file_id} EDS Measured Composition",
                template=eds_composition_measurement_template,
            )

            eds_measured_composition.update_parameters(
                Parameter(
                    "Al",
                    value=NominalReal(
                        data["Results"]["Measured Composition (%)"]["Al"], ""
                    ),
                    template=composition_parameter_template,
                ),
                Parameter(
                    "Co",
                    value=NominalReal(
                        data["Results"]["Measured Composition (%)"]["Co"], ""
                    ),
                    template=composition_parameter_template,
                ),
                Parameter(
                    "Cr",
                    value=NominalReal(
                        data["Results"]["Measured Composition (%)"]["Cr"], ""
                    ),
                    template=composition_parameter_template,
                ),
                Parameter(
                    "Cu",
                    value=NominalReal(
                        data["Results"]["Measured Composition (%)"]["Cu"], ""
                    ),
                    template=composition_parameter_template,
                ),
                Parameter(
                    "Fe",
                    value=NominalReal(
                        data["Results"]["Measured Composition (%)"]["Fe"], ""
                    ),
                    template=composition_parameter_template,
                ),
                Parameter(
                    "Mn",
                    value=NominalReal(
                        data["Results"]["Measured Composition (%)"]["Mn"], ""
                    ),
                    template=composition_parameter_template,
                ),
                Parameter(
                    "Ni",
                    value=NominalReal(
                        data["Results"]["Measured Composition (%)"]["Ni"], ""
                    ),
                    template=composition_parameter_template,
                ),
                Parameter(
                    "V",
                    value=NominalReal(
                        data["Results"]["Measured Composition (%)"]["V"], ""
                    ),
                    template=composition_parameter_template,
                ),
                which="run",
            )
            eds_measurements.append(eds_measured_composition)

            # # Create Standard Deviation Measurements for the primary phase EDS composition
            eds_measured_composition_stddev = Measurement(
                name=f"{file_id} EDS Measured Composition StdDev",
                template=eds_stddev_measurement_template,
            )

            eds_measured_composition_stddev.update_parameters(
                Parameter(
                    "Al",
                    value=NominalReal(
                        data["Results"]["Measured Composition (%)"]["Al"], ""
                    ),
                    template=composition_std_parameter_template,
                ),
                Parameter(
                    "Co",
                    value=NominalReal(
                        data["Results"]["Measured Composition (%)"]["Co"], ""
                    ),
                    template=composition_std_parameter_template,
                ),
                Parameter(
                    "Cr",
                    value=NominalReal(
                        data["Results"]["Measured Composition (%)"]["Cr"], ""
                    ),
                    template=composition_std_parameter_template,
                ),
                Parameter(
                    "Cu",
                    value=NominalReal(
                        data["Results"]["Measured Composition (%)"]["Cu"], ""
                    ),
                    template=composition_std_parameter_template,
                ),
                Parameter(
                    "Fe",
                    value=NominalReal(
                        data["Results"]["Measured Composition (%)"]["Fe"], ""
                    ),
                    template=composition_std_parameter_template,
                ),
                Parameter(
                    "Mn",
                    value=NominalReal(
                        data["Results"]["Measured Composition (%)"]["Mn"], ""
                    ),
                    template=composition_std_parameter_template,
                ),
                Parameter(
                    "Ni",
                    value=NominalReal(
                        data["Results"]["Measured Composition (%)"]["Ni"], ""
                    ),
                    template=composition_std_parameter_template,
                ),
                Parameter(
                    "V",
                    value=NominalReal(
                        data["Results"]["Measured Composition (%)"]["V"], ""
                    ),
                    template=composition_std_parameter_template,
                ),
                which="run",
            )
            eds_measurements.append(eds_measured_composition_stddev)

            # # Create Measurements for 2nd Phase EDS composition results
            eds_2nd_phase_measured_composition = Measurement(
                name=f"{file_id} 2nd Phase EDS Measured Composition",
                template=eds_composition_measurement_template,  # Reuse the same template
            )

            if "2nd Phase EDS Measured Composition (%)" in data["Results"]:
                eds_2nd_phase_measured_composition.update_parameters(
                    Parameter(
                        "Al",
                        value=NominalReal(
                            data["Results"]["2nd Phase EDS Measured Composition (%)"][
                                "Al"
                            ],
                            "",
                        ),
                        template=composition_parameter_template,
                    ),
                    Parameter(
                        "Co",
                        value=NominalReal(
                            data["Results"]["2nd Phase EDS Measured Composition (%)"][
                                "Co"
                            ],
                            "",
                        ),
                        template=composition_parameter_template,
                    ),
                    Parameter(
                        "Cr",
                        value=NominalReal(
                            data["Results"]["2nd Phase EDS Measured Composition (%)"][
                                "Cr"
                            ],
                            "",
                        ),
                        template=composition_parameter_template,
                    ),
                    Parameter(
                        "Cu",
                        value=NominalReal(
                            data["Results"]["2nd Phase EDS Measured Composition (%)"][
                                "Cu"
                            ],
                            "",
                        ),
                        template=composition_parameter_template,
                    ),
                    Parameter(
                        "Fe",
                        value=NominalReal(
                            data["Results"]["2nd Phase EDS Measured Composition (%)"][
                                "Fe"
                            ],
                            "",
                        ),
                        template=composition_parameter_template,
                    ),
                    Parameter(
                        "Mn",
                        value=NominalReal(
                            data["Results"]["2nd Phase EDS Measured Composition (%)"][
                                "Mn"
                            ],
                            "",
                        ),
                        template=composition_parameter_template,
                    ),
                    Parameter(
                        "Ni",
                        value=NominalReal(
                            data["Results"]["2nd Phase EDS Measured Composition (%)"][
                                "Ni"
                            ],
                            "",
                        ),
                        template=composition_parameter_template,
                    ),
                    Parameter(
                        "V",
                        value=NominalReal(
                            data["Results"]["2nd Phase EDS Measured Composition (%)"][
                                "V"
                            ],
                            "",
                        ),
                        template=composition_parameter_template,
                    ),
                    which="run",
                )
                eds_measurements.append(eds_2nd_phase_measured_composition)

                # # Create Standard Deviation Measurements for the 2nd Phase EDS composition
                eds_2nd_phase_measured_composition_stddev = Measurement(
                    name=f"{file_id} 2nd Phase EDS Measured Composition StdDev",
                    template=eds_stddev_measurement_template,  # Reuse the same template for standard deviation
                )
                eds_2nd_phase_measured_composition_stddev.update_parameters(
                    Parameter(
                        "Al",
                        value=NominalReal(
                            data["Results"][
                                "2nd Phase EDS Measured Composition StdDev (%)"
                            ]["Al"],
                            "",
                        ),
                        template=composition_std_parameter_template,
                    ),
                    Parameter(
                        "Co",
                        value=NominalReal(
                            data["Results"][
                                "2nd Phase EDS Measured Composition StdDev (%)"
                            ]["Co"],
                            "",
                        ),
                        template=composition_std_parameter_template,
                    ),
                    Parameter(
                        "Cr",
                        value=NominalReal(
                            data["Results"][
                                "2nd Phase EDS Measured Composition StdDev (%)"
                            ]["Cr"],
                            "",
                        ),
                        template=composition_std_parameter_template,
                    ),
                    Parameter(
                        "Cu",
                        value=NominalReal(
                            data["Results"][
                                "2nd Phase EDS Measured Composition StdDev (%)"
                            ]["Cu"],
                            "",
                        ),
                        template=composition_std_parameter_template,
                    ),
                    Parameter(
                        "Fe",
                        value=NominalReal(
                            data["Results"][
                                "2nd Phase EDS Measured Composition StdDev (%)"
                            ]["Fe"],
                            "",
                        ),
                        template=composition_std_parameter_template,
                    ),
                    Parameter(
                        "Mn",
                        value=NominalReal(
                            data["Results"][
                                "2nd Phase EDS Measured Composition StdDev (%)"
                            ]["Mn"],
                            "",
                        ),
                        template=composition_std_parameter_template,
                    ),
                    Parameter(
                        "Ni",
                        value=NominalReal(
                            data["Results"][
                                "2nd Phase EDS Measured Composition StdDev (%)"
                            ]["Ni"],
                            "",
                        ),
                        template=composition_std_parameter_template,
                    ),
                    Parameter(
                        "V",
                        value=NominalReal(
                            data["Results"][
                                "2nd Phase EDS Measured Composition StdDev (%)"
                            ]["V"],
                            "",
                        ),
                        template=composition_std_parameter_template,
                    ),
                    which="run",
                )
                eds_measurements.append(eds_2nd_phase_measured_composition_stddev)

            # Create the overall sequence for EBSD and EDS Mapping Experiment
            ebsd_eds_mapping_sequence = MaterialsSequence(
                name="EBSD and EDS Mapping Experiment",
                science_kit=science_kit,
                material=ebsd_sample_material,
                ingredients=[sample_ingredient],
                process=ebsd_eds_mapping_process,
                measurements=eds_measurements,
            )

            # Link internal elements within the sequence
            ebsd_eds_mapping_sequence.link_within()

            return ebsd_eds_mapping_sequence

        ebsd_eds_mapping_sequence = make_ebsd_eds_mapping_sequence(ebsd_eds_data)
        return science_kit

    # Function to get subfolders of a folder using the Girder API
    def get_subfolders(self, parent_folder_id):

        # API request to get subfolders of a folder
        subfolders = self.client.get(
            f"/folder",
            parameters={"parentType": "folder", "parentId": parent_folder_id},
        )

        return subfolders

    # Function to find a subfolder by its name within a parent folder
    def find_subfolder_by_name(self, parent_folder_id, folder_name):
        subfolders = self.get_subfolders(parent_folder_id)

        # Search for the folder by name
        for folder in subfolders:
            if folder["name"] == folder_name:
                return folder  # Return the folder metadata

        return None  # Return None if not found

    # Function to find a folder by traversing a given path
    def find_folder_by_path(self, root_folder_id, folder_path):
        folder_names = folder_path.strip("/").split("/")  # Split path into folder names

        current_folder_id = root_folder_id

        # Traverse through the folder hierarchy based on the folder names
        for folder_name in folder_names:
            folder = self.find_subfolder_by_name(current_folder_id, folder_name)
            if folder:
                current_folder_id = folder[
                    "_id"
                ]  # Update current folder ID to the found folder
            else:
                print(f"Folder '{folder_name}' not found.")
                return None

        return current_folder_id  # Return the ID of the final folder if found


def main(args=None):
    """
    Main method to run from command line
    """
    BIRDSHOTModeller.run_from_command_line(args)


if __name__ == "__main__":
    main()
