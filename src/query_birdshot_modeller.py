from .birdshot_modeller import BIRDSHOTModeller


class QueryBIRDSHOTModeller(BIRDSHOTModeller):

    pass


def main(args=None):
    """
    Main method to run from command line
    """
    QueryBIRDSHOTModeller.run_from_command_line(args)


if __name__ == "__main__":
    main()
