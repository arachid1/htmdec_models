import fmrest

fmrest.utils.TIMEOUT = 20
import os
import argparse
import pandas as pd
import numpy as np
import getpass


class LaserShockConstants:
    """
    Constants for Laser Shock Lab stuff and the associated FileMaker Database
    """

    @property
    def FILEMAKER_SERVER_IP_ADDRESS(self):
        return "https://10.173.38.223"  # IP Address of the FileMaker DB Server within the Hopkins VPN

    @property
    def DATABASE_NAME(self):
        return "Laser Shock"  # Name of the Laser Shock Lab's FileMaker database


LASER_SHOCK_CONST = LaserShockConstants()


def extract_metadata(fms, scope_filename):
    find_query = [{"Scope Filename": scope_filename}]
    foundset = fms.find(find_query).to_df()
    # print(foundset)
    # metadata = foundset.to_csv('{}.csv'.format(scope_filename))
    foundset.to_csv("data/filemaker_data/{}.csv".format(scope_filename))
    # metadata = foundset.to_df()
    # return metadata


def username():
    __username = None
    if __username is None:
        __username = os.path.expandvars("$JHED_UNAME")
        if __username == "$JHED_UNAME":
            __username = (input("Please enter your JHED username: ")).rstrip()
    return __username


def password():
    __password = None
    if __password is None:
        __password = os.path.expandvars("$JHED_PWORD")
        if __password == "$JHED_PWORD":
            __password = getpass.getpass(
                f"Please enter the JHED password for {username}: "
            )
    return __password


def process(scope_filenames, layout_name):
    records = {}
    fms = fmrest.Server(
        LASER_SHOCK_CONST.FILEMAKER_SERVER_IP_ADDRESS,
        user=username(),
        password=password(),
        database=LASER_SHOCK_CONST.DATABASE_NAME,
        layout=layout_name,
        verify_ssl=False,
        api_version="v1",
    )

    fms.login()

    for scope_filename in scope_filenames:
        metadata = extract_metadata(fms, scope_filename)

    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--layout_name", type=str, required=False)
    parser.add_argument("--filenames_list", type=str, required=False)
    parser.add_argument("--scope_filenames", type=str, required=False)
    args = parser.parse_args()
    # filenames_list = "data/filemaker_data/kens_velocity_flyer_scopes.txt"
    # scope_filenames = pd.read_csv(filenames_list).to_numpy()
    # scope_filenames = np.squeeze(scope_filenames)
    process(scope_filenames=args.scope_filenames, layout_name="Experiment")
    # print(fms.get_layout())
    # return records in the foundset
    # recs = fms.get_records(limit=2)
    # for rec in recs:
    #     print(rec.keys())
    #     print(rec.values())
    # final = []
    # find_query = [{'Date': '11/11/2022'}]
    # foundset = fms.find(find_query)
    # print(foundset[0].keys())
    # records_11_11 = ['00011', '00012', '00013', '00014', '00015', '00016', '000']
    # records_11_11 = list(range(11, 21)) + [26]
    # records_12_13 = list(range(5, 11)) + [21]
    # elements = [
    #     ['11/11/2022', records_11_11],
    #     ['12/13/2022', records_12_13]
    # ]
    # for i, [date, ids] in enumerate(elements):
    #     for id in ids:
    #         print(date)
    #         print(id)
    #         find_query = [{'Date':date, 'Experiment Day Counter':id}]
    #         foundset = fms.find(find_query).to_df()
    #         print(foundset)
    #         foundset.to_csv('{}.csv'.format(id))
    #         # print(foundset[0].keys())
    #         # print(foundset[0].values())
    #         if id == 12:
    #             exit()

    # for date in ['11/11/2022', '12/13/2022']:
    #     foundset = fms.find(find_query)

    # for record in foundset:
    #     print(record.values())
    # print(rec[0].__dict__)
    # print(rec.__dict__)
    # print(rec._records)


if __name__ == "__main__":
    main()
