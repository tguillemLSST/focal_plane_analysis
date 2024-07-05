import pandas as pd
import os

from lsst.summit.utils import ConsDbClient
with open("token-file", "r") as f:
    token = f.read()
client = ConsDbClient(f"https://user:{token}@usdf-rsp.slac.stanford.edu/consdb")
# Note that "user" in the URL can be any string.
print(client.schema())

##########
client.schema("lsstcomcamsim")
#client.schema("lsstcomcamsim", "cdb_lsstcomcamsim.visit1_quicklook")

instrument = 'lsstcomcamsim'
day_obs = '2024-06-25'
day_obs_int = int(day_obs.replace('-', ''))

visits_query = f'''
    SELECT * FROM cdb_{instrument}.exposure
        WHERE obs_start_mjd IS NOT NULL
        AND s_ra IS NOT NULL
        AND s_dec IS NOT NULL
        AND sky_rotation IS NOT NULL
        AND ((band IS NOT NULL) OR (physical_filter IS NOT NULL))
        AND day_obs = {day_obs_int}
'''

visits = client.query(visits_query).to_pandas()

print(visits.columns)

sys.exit()
visits

exposure_opsimdb_map = {
        'obs_start_mjd': 'observationStartMJD',
        'obs_start': 'start_date',
        's_ra': 'fieldRA',
        's_dec': 'fieldDec',
        'sky_rotation': 'rotSkyPos',
        'band': 'filter',
        'airmass': 'airmass',
        'altitude': 'altitude',
        'azimuth': 'azimuth',
        'exp_time': 'visitExposureTime'
    }

visits.rename(exposure_opsimdb_map, axis='columns', inplace=True)

missing_filter = visits['filter'].isna()
visits.loc[missing_filter, 'filter'] = visits.loc[missing_filter, 'physical_filter'].str.get(0)

if len(visits):
    displayed_columns = [
        "start_date",
        "seq_num",
        "fieldRA",
        "fieldDec",
        "filter",
        "visitExposureTime",
#        "numExposures",
#        "t_eff",
#        "skyBrightness",
#        "seeingFwhmEff",
#        "cloud",
#        "note",
    ]
    displayed_visits_df = visits.loc[:, displayed_columns]
    with pd.option_context("display.max_rows", 2000):
        display(displayed_visits_df)
else:
    print("No visits")

