import csv
import datetime
import math as m
import sys

from pathlib import Path

import archive
from formulas import vapor_pressure_liquid_water

###################################################################################################
# Get the list of files to process.
###################################################################################################
csv_dir = Path('csv')
csv_files = (x for x in csv_dir.iterdir() if x.is_file() and x.suffix == '.csv')

###################################################################################################
# Set up files for logging and KML output.
###################################################################################################
now = datetime.datetime.utcnow()
log_file = csv_dir / Path(now.strftime("%Y%m%d_%H%M%S_import.log"))
log_file = open(log_file, "wt")

kml_file = csv_dir / Path(now.strftime("%Y%m%d_%H%M%S_import.kml"))
###################################################################################################
# The variable of interest.
###################################################################################################
def vapor_pressure_deficit(t_c, td_c):
    """Given the temperature and dew point, calculate the vapor pressure deficit.

    # Arguments
    t_c - the temperature in Celsius
    td_c - the dew point in Celsius.

    # Returns
    The vapor pressure deficit in hectopascals.
    """

    saturation_vapor_pressure = vapor_pressure_liquid_water(t_c)
    vapor_pressure = vapor_pressure_liquid_water(td_c)

    return saturation_vapor_pressure - vapor_pressure

###################################################################################################
# Map the QC codes in the data provided NCEI to my own, more mathematical system.
###################################################################################################
#
# Since the variable I'm interested in (Vapor Pressure Deficit, VPD) is a function of temperature
# and dew point, I'll use the worst case for each of those codes when determining the code for my
# VPD values.
#
# Lower values are better, so if I have two measurements for a valid time and station, I'll want to
# keep the one with the lower value for the QC code.
#
# I'll also use a global dictionary to keep track of how often each of the QC codes from the source
# data shows up so I have a good understanding of how the data has been processed before I got it.
#
# From the documentation I'll list what each code means here:
# 0 = Passed gross limits check
# 1 = Passed all quality control checks
# 2 = Suspect
# 3 = Erroneous
# 4 = Passed gross limits check, data originate from an NCEI data source
# 5 = Passed all quality control checks, data originate from an NCEI data source
# 6 = Suspect, data originate from an NCEI data source
# 7 = Erroneous, data originate from an NCEI data source
# 9 = Passed gross limits check if element is present
# A = Data value flagged as suspect, but accepted as a good value
# C = Temperature and dew point received from Automated Weather Observing System (AWOS) are reported 
#     in whole degrees Celsius. Automated QC flags these values, but they are accepted as valid.
# I = Data value not originally in data, but inserted by validator
# M = Manual changes made to value based on information provided by NWS or FAA
# P = Data value not originally flagged as suspect, but replaced by validator
# R = Data value replaced with value computed by NCEI software
# U = Data value replaced with edited value
qc_codes_counts = {
        '0': 0, '1': 0, '2': 0, '3': 0, '4': 0, '5': 0, '6': 0, '7': 0, '9': 0,
        'A': 0, 'C': 0, 'I': 0, 'M': 0, 'P': 0, 'R': 0, 'U': 0, 
}


def map_qc_codes(tmp_qc, dew_qc):
    '''Take the temperature and dew point QC codes and combine them into a single code.

    Lower is better.

    # Returns
    an integer.
    '''
    def map_single_code(qc):
        # Update the global dictionary of counts.
        qc_codes_counts[qc] += 1
        if qc == '5':
            return 0
        elif qc in set(['0', '1', '4', '9', 'A', 'I', 'M', 'P', 'R', 'U']):
            return 5
        elif qc == 'C':
            return 10
        elif qc in set(['2', '6']):
            return 20
        elif qc in set(['3', '7']):
            return 50
        else:
            raise Exception(f"Unknown QC code {qc}")

        return

    return max(map_single_code(tmp_qc), map_single_code(dew_qc))


def output_qc_codes_counts():
    '''Write the content of qc_codes_counts to the global log file.'''

    log_file.write("\nNumber of occurrences of QC codes.\n")
    for code, count in qc_codes_counts.items():
        log_file.write(f"{code:4s} {count:10d}\n")

###################################################################################################
# Convert rows of the input data (csv files) into a consistent format.
###################################################################################################
def build_row_converter(cols):
    '''Build a function to convert rows from a CSV file into a common format.

    # Arguments
    cols - is the first line of the CSV file and contains the column names.

    # Returns
    A function specifically tailored to parsing a file with the provided columns in the provided 
    order. The returned function takes a row from the CSV module and parses the columns into the
    desired format. In the event that a certain column can't be parsed, it will return `None`.
    '''

    # Find the column indexes of the needed values.
    stn_idx = cols.index('STATION')

    name_idx = cols.index('NAME')
    lat_idx = cols.index('LATITUDE')
    lon_idx = cols.index('LONGITUDE')
    elev_idx = cols.index('ELEVATION')

    id_idx = cols.index('CALL_SIGN')
    date_idx = cols.index('DATE')
    tmp_idx = cols.index('TMP')
    dew_idx = cols.index('DEW')

    def converter(row):
        try:
            stn_num = row[stn_idx].strip()
            name = row[name_idx].strip()

            lat = float(row[lat_idx].strip())
            lon = float(row[lon_idx].strip())
            elev = float(row[elev_idx].strip())
            coords = (lat, lon, elev)

            stn_id = row[id_idx].strip().lower()
            if '9999' in stn_id or stn_id == '':
                stn_id = None
            vt = row[date_idx].strip()

            tmp = row[tmp_idx].strip()
            dew = row[dew_idx].strip()

            if '9999' in tmp or '9999' in dew:
                # These are the missing value indicators for the temperature and dew point. If we're
                # missing either of these variables on a row, nothing on the row is useful to us.
                return None

            stn_num = int(stn_num)

            vt = datetime.datetime.strptime(vt, "%Y-%m-%dT%H:%M:%S")

            tmp, tmp_qc = tmp.split(',')
            dew, dew_qc = dew.split(',')

            tmp = float(tmp) / 10
            dew = float(dew) / 10
            vpd = vapor_pressure_deficit(tmp, dew)
            qc_code = map_qc_codes(tmp_qc, dew_qc)

            return (stn_num, name, stn_id, coords, vt, tmp, dew, vpd, qc_code)

        except ValueError:
            # For any parsing error, just skip the whole row.
            return None

    return converter

###################################################################################################
# Round a time to the nearest hour.
###################################################################################################
def round_to_hour(vt):
    hour = vt.hour
    minute = vt.minute
    year, month, day, hour, *_ = vt.timetuple()

    if minute > 30:
        rounded = vt + datetime.timedelta(60 - minute)
    else:
        rounded = vt - datetime.timedelta(minute)

    year, month, day, hour, *_ = rounded.timetuple()
    rounded = datetime.datetime(year, month, day, hour)

    return rounded

###################################################################################################
# Stations to skip and not import into the database.
###################################################################################################
#
# Most stations are listed to skip because they are missing too much data to be useful for long term
# analysis. This isn't needed because I (should) have data checks at the time I query the data, but
# I've left it because it keeps the database smaller, and consequently helps programs that query it
# run faster.
#
# I add stations to this list that have so little data they really have no chance of being used in 
# any analysis I have in mind. It would be a good idea to reevaluate the stations in this this list 
# if I reuse this code for another project.
skip_stations = {
        # Canada
        71783099999: 'cwyj', # VICTORIA UNIVERSITY CS, CA
        71799499999: 'cwyj',

        71914099999: 'cvts', # SATURNA CAPMON CS, CA

        71798099999: 'cwpf', # ESQUIMALT HARBOUR

        71780099999: 'cwsp', # kwsp, SHERINGHAM POINT BC, CA
        71200499999: 'cwsp',

        71778099999: 'cwqk', # RACE ROCKS CAMPBELL SCIENTIFIC BC, CA

        71473599999: 'cywh', # VICTORIA HARBOUR, CA

        71473099999: 'kwez', # SATURNA ISLAND CS, CA

        71799599999: 'kwir', # VICTORIA MARINE RAD
        71202399999: 'kwir', 
        72799599999: 'kwir',

        71483399999: 'kwkh', # MALAHAT AUTO8
        71774099999: 'kwkh', 

        71200399999: 'kwlm', # VICTORIA AUTO8
        71200099999: 'kwlm',

        71202699999: 'kwpf', # ESQUIMALT METOC

        74200099999: 'kwvg', # VICTORIA AUTO8

        71202799999: 'kwvv', # VICTORIA HARTLAND CS

        71031099999: 'kwzv', # DISCOVERY ISLAND

        71927099999: 'kyur', # NORTH COWICHAN

        71034099999: 'xtic', # TRIAL ISLAND

        # Washington
        74200699999: 'k75s', # BURLINGTON MOUNT VERN

        72796499999: 'k76s', # (kokh) OAK HARBOR

        72793599999: 'kbif', # Boeing Field

        72027299999: 'kbvs', # SKAGIT REGIONAL
        99999994282: 'kbvs',
        72027294282: 'kbvs',

        72025400119: 'kcls', # CHEHALIS CENTRALIA AIRPORT, WA US
        94999900119: 'kcls',

        72787399999: 'kcqv', # COLVILLE MUNICIPAL AIRPORT
        72787199999: 'kcqv',

        72781399999: 'kfct', # Vagabond Army Air Field, Yakima

        74207099999: 'kgrf', # Fort Lewis near Gray AAF

        72784094187: 'khms', # Hanford
        72784099999: 'khms',
        99999994187: 'khms',

        74201024228: 'knou', # 'know', 'knkw' PORT ANGELES WEATHER BUREAU AIRPORT
        99999924228: 'knou',

        72789099999: 'komk', # Omak

        72220804224: 'kors', # EASTSOUND ORCAS ISLAND AIRPORT
        99999904224: 'kors',
        72220899999: 'kors',

        72038800469: 'kplu', # PIERCE CO AIRPORT THUN FIELD

        72788499999: 'krld', # Richland Airport

        99999994248: 'krnt', # RENTON MUNICIPAL

        74201699999: 'ks19', # FRIDAY HARBOR

        72796599999: 'ks88', # ARLINGTON AWOS
        74918599999: 'ks88',

        72792699999: 'ktdo', # ED CARLSON MEM FIELD
        72792624241: 'ktdo',
        99999924241: 'ktdo',

        99430099999: 'ktti', # TATOOSH ISLAND WA
        99999924240: 'ktti',

        99999924252: 'kuil', # QUILLAYUTE 

        99407099999: 'xdsw', # DESTRUCTION IS. WA

        99999924278: 'xhoq', # HOQUIAM

        99999924147: 'xlac', # La Crosse

        99999924226: 'xnhw', # NORTH HEAD

        74201099999: 'xpaw', # PORT ANGELES CGAS

        99800799999: 'xpbf', # PADILLA BAY RESERVE, WA US

        99418099999: 'xsis', # SMITH ISLAND WA

        99999924253: 'xsnw', # SHELTON NAAS

        99999924238: 'xstv', # STEVENSON

        99435099999: 'xwpw', # WEST POINT WA, US

        72056299999: 'xyak', # RANGE OP 13 YAKIMA TRAINING CENTER

        # Oregon
        99424099999: 'caro', # CAPE ARAGO OR

        99206099999: 'db20', # ENVIRONM BUOY 46041

        99237099999: 'db23', # COL RIVER BAR 78NM SOUTH SOUTHWEST OF

        99246099999: 'db24', # DATA BUOY 46010

        99275099999: 'db27', # DATA BUOY 46043

        99297099999: 'db29', # ENVIRONM BOUY 46039

        99631099999: 'db63', # MOORED BUOY 46046

        72589699999: 'k4lw', # LAKEVIEW
        72589199999: 'k4lw',
        72589099999: 'k4lw',

        72683994196: 'k5j0', # JOHN DAY AIRPORT
        99999994196: 'k5j0',
        72687600387: 'k5j0',
        72687699999: 'k5j0',
        72683999999: 'k5j0',

        72691699999: 'k90s', # UMPQUA RIVER COAST GUARD STATION

        72063800224: 'kbdn', # BEND MUNICIPAL AIRPORT
        94999900224: 'kbdn',

        72598524267: 'kbok', # BROOKINGS STATE AIRPORT
        72598599999: 'kbok',
        72598199999: 'kbok',
        72598099999: 'kbok',
        72036524267: 'kbok',
        72036599999: 'kbok',

        99999924134: 'kbno', # Federal Building, Burns, OR

        72694599999: 'kcvo', # CORVALLIS MUNICIPAL
        99999924248: 'kcvo',

        72698799999: 'kczk', # CASCADE LOCKS STATE
        99999994204: 'kczk',
        72698794204: 'kczk',

        99999924249: 'klmt', # Klamath

        72695199999: 'konp', # Newport
        99428099999: 'konp',

        72687499999: 'kp88',  # ROME AUTOMATIC METEOROLOGICAL OB

        72695799999: 'kpfc', # PACIFIC CITY STATE

        74918499999: 'ks33', # CITY CO AIRPORT

        72683604201: 'kspb', # SCAPPOOSE INDUSTRIAL AIRPORT
        99999904201: 'kspb',
        72683699999: 'kspb',

        72020200118: 'ktmk', # TILLAMOOK AIRPORT
        72696399999: 'ktmk',
        94999900118: 'ktmk',
        72020299999: 'ktmk',
        99999924254: 'ktmk',

        72599099999: 'xblnc', # CAPE BLANCO

        99799599999: 'xchn', # SOUTH SLOUGH RESERVE

        72083799999: 'xjos', # JOSEPH STATE

        99999994285: 'xklkv', # LAKEVIEW LAKE CO AIRPORT? Large variation in coords.

        99999924250: 'xlvn', # LAKEVIEW NSA

        99999924251: 'xnbn', # NORTH BEND NAAF

        99999924236: 'xsky', # SISKIYOU SUMMIT

        72688799999: 'xsrv', # SUNRIVER, OR US

        # Idaho
        72686099999: 'k27u', # Salmon, but 46 km away?
        72686699999: 'k27u',

        72578724158: 'k4sv', # STREVELL
        99999924158: 'k4sv',

        72092300310: 'k65s', # Bonner's Ferry
        94999900310: 'k65s',

        72586999999: 'k77m', # MALTA
        72586994186: 'k77m',

        72680999999: 'kcbx', # BOISE NEXRAD

        99999924140: 'kdbs', # DUBOIS FAA AIRPORT

        72036904135: 'kdij', # DRIGGS REED MEMORIAL AIRPORT

        72687300452: 'kgic', # Grangeville
        72687399999: 'kgic',
        72687199999: 'kgic',

        72586099999: 'kgng', # GOODING MUNICIPAL
        72586124142: 'kgng',

        72073400264: 'kman', # NAMPA MUNICIPAL AIRPORT
        94999900264: 'kman',

        72578624151: 'kmld', # MALAD CITY
        72578699999: 'kmld',
        99999924151: 'kmld',

        72783699999: 'kmlp', # Mullan Pass, this location is down in town.

        72587199999: 'kowy', # Middle of nowhere. {42.583} {-116.167} {', US'}

        72450399999: 'ks14', # SPENCER

        72687099999: 'ks80', # Grangeville, but about 3.6 km away.

        72682404112: 'ksnt', # STANLEY RANGER STATION

        99999924196: 'ksmn', # Salmon

        72586599999: 'ksun', # FRIEDMAN MEM
        99999994161: 'ksun',
        72586594161: 'ksun',

        72032299999: 'kszt', # Sandpoint
        99999904129: 'kszt',
        72032204129: 'kszt',

        72586894184: 'ku78', # ALLEN H TIGERT AIRPORT
        72586899999: 'ku78',
        99999994184: 'ku78',

        99999924142: 'xgoo', # GOODING 2 S

        99999994143: 'xidf', # IDAHO FALLS 46 W

        74413099999: 'xsay', # SAYLOR CREEK GUNNERY RANGE ID.

        # Montana
        72099599999: 'k1bm', #  BRAVO GEYSER

        72099699999: 'k1cm', #  CHARLIE STANFORD

        72099999999: 'k1fm', #  FOXTROT AUGUSTA

        72100399999: 'k1im', # India ULM

        72100499999: 'k1jm', #  JULIE POWER

        72100599999: 'k1km', # Kilo Harlowtown

        72100899999: 'k1nm', # NOVEMBER GRASS RANGE

        72773524139: 'k3du', # Drummond
        72773599999: 'k3du',
        99999924139: 'k3du',

        72677599999: 'k3ht', # Harlowtown

        72779599999: 'k3th', #  THOMPSON FALLS

        72667099999: 'k4bq', # Broadus
        72667199999: 'k4bq',
        72667399999: 'k4bq',

        72679524161: 'k4ha', # Whitehall
        99999924161: 'k4ha',

        72075399999: 'k4s2', # KEN JERNSTEDT

        74918399999: 'k79s', # FORT BENTON

        72679099999: 'kbtm2', # Butte, but 2.5km away.

        72099899999: 'kem1', # Echo Winifred
        
        72772599999: 'ketm', # Elliston

        99999994010: 'kgsg', #  ST. MARIE

        72768599999: 'kgsg2', #  GLASGOW INDUSTRIAL

        72776099999: 'kgtf', # Great Falls

        72677499999: 'khmm', #  HAMILTON RAVALLI CO
        72075999999: 'khmm',

        99999924035: 'khvr', # Havre Weather Office

        72100699999: 'klm1', #  LIMA JUDITH GAP

        72100799999: 'kmm1', #  MIKE MOORE

        72676699999: 'kmqm', # Monida

        74918299999: 'kp93', #  JUDITH GAP

        72676494163: 'kwys', # West Yellowstone
        99999994163: 'kwys',
        72676499999: 'kwys',
        72676599999: 'kwys',

        72676399999: 'kwys2', # West Yellowstone, but different.
        99999924198: 'kwys2',
        72676324198: 'kwys2',
        72676199999: 'kwys2',

        72099499999: 'xarn', #  ALPHA RAYNESFORD

        99999924040: 'xcus', #  CUSTER

        72099799999: 'xded', #  DELTA DENTON

        99999924034: 'xggw', #  GLASGOW NUMBER 2

        72100199999: 'xgsm', #  GOLF SIMMS

        72100299999: 'xhff', #  HOTEL FAIRFIELD

        72100999999: 'xosr', #  OSCAR ROY

        94999900114: 'xsba', #  STARR BROWNING AIRSTRIP

        99999924159: 'xsup', #  SUPERIOR

        99999994045: 'xxmt', #  SOMEWHERE IN SE MT

        # Wyoming
        72568524019: 'k4dg', # DOUGLAS
        72568099999: 'k4dg',
        72568124019: 'k4dg',
        99999924019: 'k4dg',

        72663099999: 'k4mc', # MOORCROFT
        99999924088: 'k4mc', 

        72577199999: 'kbpi2', # BIG PINEY - but about 5km away?
        72577099999: 'kbpi2', # This location moved about 5km without being renamed.
        72671099999: 'kbpi',

        99999924045: 'kcod2', # Cody but about 12km away, big range in coords for this location.

        72568699999: 'kdgw', # CONVERSE CO
        72568694057: 'kdgw',
        99999994057: 'kdgw',

        72052200443: 'kdub', # DUBOIS MUNICIPAL AIRPORT

        99999900446: 'kecs', # MONDELL FIELD AIRPORT

        72097199999: 'kfwz', # SOUTH PASS
        72097100342: 'kfwz',
        94999900342: 'kfwz',

        72096900341: 'kjpd', # TEN SLEEP
        72096999999: 'kjpd',
        94999900341: 'kjpd',

        72034599999: 'kpna', # RALPH WENZ FIELD
        72034594086: 'kpna',

        72564899999: 'ksib', # SIBLEY PEAK

        72103600353: 'ktbx', # BOYSEN THERMOPOL
        72103699999: 'ktbx',
        94999900353: 'ktbx',

        72576394053: 'ktor', # TORRINGTON MUNICIPAL AIRPORT
        99999994053: 'ktor',
        72576399999: 'ktor',

        99999900484: 'kw43', # HULETT MUNICIPAL AIRPORT

        72676099999: 'kwey', # West Yellowstone, but different, big range in coords

        72051699999: 'xafm', # AFTON MUNICIPAL

        99999924016: 'xcas', # CASPER

        72068499999: 'xphi', # PHIFER AFLD

        72051999999: 'xpwl', # POWELL MUNICIPAL

        # North Dakota
        72101899999: 'k1hn', # Hotel Parshall

        72101999999: 'k1in', # India Palermo

        72102299999: 'k1kn', # KILO DONNYBROOK 2

        72102399999: 'k1ln', # LIMA BOWBELLS

        72529699999: 'kbpp', # Bowman

        72086100287: 'kd50', # Crosby

        99999924014: 'kisn', # Weather Office

        72086899999: 'ks25', # Watford 
        72086800294: 'ks25',
        94999900294: 'ks25',

        99999994099: 'kxwa', # WILLISTON

        94999900287: 'xcma', # CROSBY MUNICIPAL AIRPORT

        # South Dakota
        72090900300: 'k08d', # Stanley
        94999900300: 'k08d', 

        72085400282: 'k20u', # Beach
        94999900282: 'k20u',

        72086300289: 'kd60', # Tioga
        94999900289: 'kd60',

        72660699999: 'kefc', # Belle Fourche Municiple

        72661099999: 'krej', # REDIG SD.

        72660500386: 'kspf', # BLACK HILLS AIRPORT CLYDE ICE FIELD
        72660599999: 'kspf',

        72669199999: 'ky22', # Lemmon
        72669599999: 'ky22',

        # Nebraska
        72563799999: 'kgrn', # GORDON MUNICIPAL

}

###################################################################################################
# Mapping of station numbers to staion ids.
###################################################################################################
#
# For various reasons (that I do not fully understand), the data from a specific site may have 
# different station numbers, usually covering different date ranges. This map takes different 
# station numbers and maps them to a list of sites, so that data from different station numbers but
# the same geographical location can be grouped together.
#
# The site id's I used in here are arbitrary. I tried to use values already associated with the 
# location such as any ICAO identifiers listed in the 'CALL SIGN' column of the source data. But in
# some cases I just had to make something up. The ones I made up usually start with 'x'.
stations = {
        # Canada - yes, a few made it in.
        71799099999: 'kyyj', # VICTORIA INTERNATIONAL
        72799099999: 'kyyj', 

        # Washington
        72784624160: 'kalw', # Walla Walla
        99999924160: 'kalw',

        72794599999: 'kawo', # ARLINGTON MUNICIPAL
        99999904205: 'kawo',
        72794504205: 'kawo',

        72793524234: 'kbif', # 'kbfi' SEATTLE BOEING FIELD
        99999924234: 'kbif',

        72797624217: 'kbli', # BELLINGHAM INTERNATIONAL AIRPORT

        72788594266: 'kclm', # PORT ANGELES FAIRCHILD INTERNATIONAL AIRPORT, WA US
        72788599999: 'kclm',
        99999994266: 'kclm',

        99999924207: 'ktcm', # TACOMA MCCHORD AFB

        72785499999: 'kdew', # Deer Park
        72787099999: 'kdew',
        99999994119: 'kdew',
        72787094119: 'kdew',
        72785494119: 'kdew',

        72698824219: 'kdls', # The Dalles
        99999924219: 'kdls',

        72782599999: 'keat', # Pangborn Memorial
        72782594239: 'keat',

        72788399999: 'keln', # Ellensburg Bowers Field
        72788324220: 'keln',
        99999924220: 'keln',

        72782624141: 'keph', # Ephrata
        72790024141: 'keph',
        99999924141: 'keph',

        72798594276: 'kfhr', # FRIDAY HARBOR AIRPORT, WA US
        72798599999: 'kfhr',
        99999994276: 'kfhr',

        72785024157: 'kgeg', # Spokane International, Geiger Field
        99999924157: 'kgeg',

        74207024201: 'kgrf', # GRAY ARMY AIR FIELD
        74207124201: 'kgrf',
        74207199999: 'kgrf',

        72792394225: 'khqm', # HOQUIAM BOWERMAN AIRPORT, WA US
        72792794225: 'khqm',
        99999994225: 'khqm',

        72792424223: 'kkls', # KELSO SOUTHWEST REGIONAL AIRPORT
        99999924223: 'kkls',
        72792499999: 'kkls',

        72782799999: 'kmwh', # Moses Lake
        72782724110: 'kmwh',

        99999924244: 'knej', # SEATTLE NAS

        72797524255: 'knuw', # WHIDBEY ISLAND NAS
        69023024255: 'knuw',
        72074924255: 'knuw',
        99999924255: 'knuw',

        72792024227: 'kolm', # OLYMPIA AIRPORT, WA US
        99999924227: 'kolm',

        72789094197: 'komk', # Omak
        99999994197: 'komk',

        72793724222: 'kpae', # EVERETT SNOHOMISH CO AIRPORT
        72793799999: 'kpae',
        99999924222: 'kpae',

        72784524163: 'kpsc', # Pasco Tri Cities
        72784599999: 'kpsc',
        99999924163: 'kpsc',

        72785799999: 'kpuw', # Pullman (WA) / Moscow (ID)
        99999994129: 'kpuw',
        72785794129: 'kpuw',

        72792894263: 'kpwt', # BREMERTON, WA US
        72792899999: 'kpwt', 
        99999994263: 'kpwt',

        72793499999: 'krnt', # RENTON MUNICIPAL
        72793494248: 'krnt',

        72793024233: 'ksea', # SEATTLE TACOMA AIRPORT, WA US
        99999924233: 'ksea',

        72785699999: 'ksff', # Felts Field
        99999994176: 'ksff',
        72785694176: 'ksff',
        72785694176: 'ksff',

        72792594227: 'kshn', # SHELTON AIRPORT
        72792599999: 'kshn',

        72785524114: 'kska', # Fairchild AFB
        72785599999: 'kska',
        69148499999: 'kska',

        72781524237: 'ksmp', # Stampede Pass
        99999924237: 'ksmp',

        74206024207: 'ktcm', # TACOMA MCCHORD AFB, WA US
        74206099999: 'ktcm',

        72793894274: 'ktiw', # TACOMA NARROWS AIRPORT, WA US
        72793899999: 'ktiw',
        99999994274: 'ktiw',

        72797094240: 'kuil', # QUILLAYUTE AIRPORT
        99999994240: 'kuil',

        72781024243: 'kykm', # Yakima Airport
        99999924243: 'kykm',

        72791894298: 'kvuo', # VANCOUVER PEARSON AIRPORT, WA US
        72791899999: 'kvuo',
        99999994298: 'kvuo',

        72698599999: 'kttd', # PORTLAND TROUTDALE
        99999924242: 'kttd',
        72698524242: 'kttd',

        # Oregon
        72791094224: 'kast', # ASTORIA AIRPORT PORT OF
        99999994224: 'kast',
        99999924246: 'kast',

        72688624130: 'kbke', # Baker City
        99999924130: 'kbke',

        72683094185: 'kbno', # BURNS MUNICIPAL AIRPORT
        72683894185: 'kbno',
        72683099999: 'kbno',

        72694524202: 'kcvo', # CORVALLIS MUNICIPAL
        99999924202: 'kcvo',

        72693024221: 'keug',  # EUGENE MAHLON SWEET AIRPORT
        99999924221: 'keug',

        72698694261: 'khio', # PORTLAND HILLSBORO AIRPORT, OR US
        72698699999: 'khio',
        99999994261: 'khio',

        72688399999: 'khri', # Hermiston 
        72688304113: 'khri',
        99999904113: 'khri',

        72688499999: 'klgd', # La Grande Union
        72688424148: 'klgd',
        99999924148: 'klgd',

        72597694285: 'klkv', # LAKEVIEW LAKE CO AIRPORT
        72597699999: 'klkv',

        72589599999: 'klmt', # KLAMATH FALLS
        72589594236: 'klmt',
        99999994236: 'klmt',
        99999924224: 'klmt',

        72688524152: 'kmeh', # Meacham
        99999924152: 'kmeh',

        72597024225: 'kmfr', # MEDFORD INTERNATIONAL AIRPORT
        99999924225: 'kmfr',

        72688199999: 'kmmv', # MCMINNVILLE MUNICIPAL
        99999994273: 'kmmv',
        72688194273: 'kmmv',

        72695024285: 'konp', # NEWPORT MUNICIPAL AIRPORT
        72695899999: 'konp',
        72695499999: 'konp',
        72695824285: 'konp',
        72695099999: 'konp',
        99999924285: 'konp',

        72691799999: 'koth', # NORTH BEND MUNICIPAL
        72691724284: 'koth',
        72691199999: 'koth',
        99999924284: 'koth',
        72691099999: 'koth',

        72688024155: 'kpdt', # Pendleton
        99999924155: 'kpdt',

        72698024229: 'kpdx', # PORTLAND INTERNATIONAL AIRPORT
        99999924229: 'kpdx',

        72690424231: 'krbg', # ROSEBURG REGIONAL AIRPORT
        72690499999: 'krbg',
        72690124231: 'krbg',
        99999924231: 'krbg',

        72692024230: 'krdm', # REDMOND AIRPORT
        72683524230: 'krdm',
        99999924230: 'krdm',

        72687599999: 'kreo', # ROME STATE
        72687594107: 'kreo',
        99999994107: 'kreo',

        72694024232: 'ksle', # SALEM AIRPORT MCNARY FIELD, OR US
        99999924232: 'ksle',

        72695994281: 'kuao', # AURORA STATE AIRPORT
        72695999999: 'kuao',
        99999994281: 'kuao',

        # Idaho
        72681024131: 'kboi',  # BOISE AIR TERMINAL
        99999924131: 'kboi',

        72586724133: 'kbyi', # BURLEY MUNICIPAL AIRPORT
        99999924133: 'kbyi',

        72783499999: 'kcoe', # Coure D'Alene
        72783424136: 'kcoe',
        99999924136: 'kcoe',

        72681399999: 'keul', # CALDWELL AWOS
        72681394195: 'keul', 
        99999994195: 'keul',

        72578524145: 'kida', # IDAHO FALLS FANNING FIELD
        99999924145: 'kida',

        72681604110: 'kjer', # JEROME CO AIRPORT
        99999904110: 'kjer',
        72681699999: 'kjer',

        72214204114: 'kllj', # Challis
        72214299999: 'kllj',
        72783399999: 'kllj',
        72783304114: 'kllj',

        72783024149: 'klws', # Lewiston
        99999924149: 'klws',

        72681724154: 'kmlp', # Mullan Pass
        72681799999: 'kmlp', 
        99999924154: 'kmlp', 

        72681599999: 'kmuo', # MOUNTAIN HOME AFB
        72681524106: 'kmuo', 

        99999994182: 'kmyl', # McCall
        72578899999: 'kmyl',
        72578894182: 'kmyl',
        72586499999: 'kmyl',
        72586494182: 'kmyl',

        72578499999: 'kp69', # Powell
        99999904109: 'kp69',
        72578404109: 'kp69',

        72578024156: 'kpih', # POCATELLO REGIONAL AIRPORT
        99999924156: 'kpih',

        72681894194: 'krxe', # REXBURG MADISON CO AIRPORT

        72686599999: 'ksmn', # Salmon
        72686199999: 'ksmn',
        72686524196: 'ksmn',

        72586694178: 'ktwf', # TWIN FALLS SUN VALLEY REGIONAL AIRPORT
        99999994178: 'ktwf',
        72586699999: 'ktwf', 

        # Montana
        72677794055: 'kbhk', # Baker
        72677799999: 'kbhk',
        99999994055: 'kbhk',

        72677024033: 'kbil', # Billings
        99999924033: 'kbil',

        99999924135: 'kbtm', # Butte
        72774024135: 'kbtm',
        72678524135: 'kbtm',
        72679124135: 'kbtm',

        72679724132: 'kbzn', # Bozeman
        99999924132: 'kbzn',

        99999924137: 'kctb', # Cut Bank
        72779624137: 'kctb',
        72769024137: 'kctb',
        
        72770024138: 'kdln', # Dillon
        72679624138: 'kdln',
        72679999999: 'kdln',
        99999924138: 'kdln',

        72667624087: 'kgdv', # Glendive
        99999924087: 'kgdv',
        72667699999: 'kgdv',

        72775524112: 'kgfa', # Malstrom AFB
        72775599999: 'kgfa',
        
        72768094008: 'kggw', # Glasgow
        99999994008: 'kggw',

        99999924146: 'kgpi', # Kalispell
        72779024146: 'kgpi',
        
        99999924143: 'kgtf', # Great Falls
        72775024143: 'kgtf',
        
        99999924144: 'khln', # Helena
        72772024144: 'khln',

        72777094012: 'khvr', # Havre
        99999994012: 'khvr',

        99999994051: 'kjdn', # Jordan
        72768494051: 'kjdn',
        72768499999: 'kjdn',
        
        72677624036: 'klwt', # Lewistown
        99999924036: 'klwt',

        72679824150: 'klvm', # Livingston
        99999924150: 'klvm',

        74230024037: 'kmls', # Miles City
        72667524037: 'kmls',
        99999924037: 'kmls',

        72773024153: 'kmso', # Missoula
        99999924153: 'kmso',

        72768694017: 'kolf', # Wolf Point
        99999994017: 'kolf',
        72768699999: 'kolf',

        72768794028: 'ksdy', # Sidney
        72768799999: 'ksdy',
        99999994028: 'ksdy',

        # Wyoming
        72671024164: 'kbpi', # BIG PINEY MARBLETON AIRPORT

        72665494054: 'kbyg', # BUFFALO JOHNSON CO AIRPORT
        99999994054: 'kbyg',
        72665499999: 'kbyg',
        72670499999: 'kbyg',

        72670099999: 'kcod', # Cody
        72670024045: 'kcod',
        72666699999: 'kcod',

        72569024089: 'kcpr', # CASPER NATRONA CO INTERNATIONAL AIRPORT
        99999924089: 'kcpr',

        72665099999: 'kgcc', # GILLETTE GILLETTE C
        99999994023: 'kgcc',
        72665094023: 'kgcc',
        72666299999: 'kgcc',
        72665599999: 'kgcc',

        72666799999: 'kgey', # South Big Horn Co
        99999924048: 'kgey',
        72666724048: 'kgey',

        72577624166: 'kjac', # JACKSON HOLE AIRPORT
        99999924166: 'kjac',
        72577699999: 'kjac',

        72576024021: 'klnd', # LANDER HUNT FIELD AIRPORT
        99999924021: 'klnd',

        72666494173: 'kp60', # Yellowstone Lake

        72576524061: 'kriw', # RIVERTON REGIONAL AIRPORT
        72672099999: 'kriw',
        72576599999: 'kriw',
        72672024061: 'kriw',

        72666024029: 'kshr', # Sheridan
        99999924029: 'kshr',

        72666524062: 'kwrl', # WORLAND MUNICIPAL AIRPORT
        72666599999: 'kwrl',

        # North Dakota
        72764524012: 'kdik', # Dickinson
        72763024012: 'kdik',
        99999924012: 'kdik',

        72758494038: 'khei', # Hettinger
        72758499999: 'khei',
        99999994038: 'khei',

        72767094014: 'kisn', # Williston
        99999994014: 'kisn', 

        # South Dakota
        72662794037: 'k2wx', # Buffalo
        72662799999: 'k2wx',

        72651499999: 'kcut', # CUSTER CO
        72651494032: 'kcut',

        99999994056: 'kd07', # Faith Municipal
        72653999999: 'kd07', 
        72653994056: 'kd07',

        72651794039: 'kien', # PINE RIDGE AIRPORT
        72651799999: 'kien',
        99999994039: 'kien',

        72683724162: 'kono', # ONTARIO
        72683799999: 'kono',
        99999924162: 'kono',

        72662024090: 'krap', # RAPID CITY REGIONAL AIRPORT
        99999924090: 'krap',

        72662599999: 'krca', # ELLSWORTH AFB
        72662524006: 'krca',
        99999924026: 'krca',

        72597524235: 'ksxt', # SEXTON SUMMIT
        99999924235: 'ksxt',

        # Nebraska !
        72563599999: 'kaia', # ALLIANCE MUNICIPAL
        72563524044: 'kaia',
        99999924044: 'kaia',

        72563624017: 'kcdr', # CHADRON MUNICIPAL AIRPORT
        72563699999: 'kcdr',
        99999924017: 'kcdr',

}

###################################################################################################
# Import the data from a file.
###################################################################################################
def import_data_from_csv(fname):
    '''Import data from a given file.
    
    This function can be called multiple times and each time it will read from the global 
    dictionary `skip_stations` and modify the `stations` dictionary. 

    # Arguments
    fname - is the name or path to the file to process.

    # Returns
    This function creates a generator object that yeilds tuples of (station number, station id, 
    coordinates as (lat, lon, elevation), valid time, VPD, QC code value).
    '''
    with open(fname, 'r') as f:
        csvreader = csv.reader(f, delimiter=',', quotechar='"')

        cols = next(csvreader)
        
        convert_row = build_row_converter(cols)

        data = (convert_row(r) for r in csvreader)
        data = (r for r in data if r is not None)

        for stn_num, name, stn_id, coords, vt, tmp, dew, vpd, qc_code in data:

            if stn_num in stations:
                stn_id_forced = stations[stn_num]
                if stn_id_forced != stn_id:
                    name = f"{stn_id} {name}"
                    stn_id = stn_id_forced

            vt = round_to_hour(vt)

            yield (stn_num, stn_id, coords, vt, tmp, dew, vpd, qc_code, name)

    return


###################################################################################################
# The main part of the script that puts it all together.
###################################################################################################

def main():
    do_prune = '-p' in sys.argv
    do_averages = '-a' in sys.argv
    do_import = "-i" in sys.argv

    log_file.write(f"  prune = {do_prune}\naverage = {do_averages}\n  import = {do_import}\n")
    print(f"  prune = {do_prune}\naverage = {do_averages}\n  import = {do_import}")

    if not (do_prune or do_averages or do_import):
        print("No operation selected. Choose -p and/or -a and/or -i")

    # Open an archive and process each file.
    with archive.Archive() as ar:
        # Make sure and remember which stations to skip.
        ar.add_skip_stations(skip_stations.keys())

        # Add new data
        if do_import:
            for fname in csv_files:
                log_file.write(f"Loading {fname}\n")
                print(f"Loading {fname}")

                data_gen = import_data_from_csv(fname)
                ar.add_all(data_gen)

        # Remove data we don't want.
        if do_prune:
            num_removed = ar.remove_stations(skip_stations.keys())
            log_file.write(f"Removed {num_removed} observations.\n")
            print(f"Removed {num_removed} observations.")

        # pre-build multi-hour averages.
        if do_averages:
            ar.build_xhour_vpd_averages()

        ar.output_sites_kml(kml_file, log_file)

    if do_import:
        # Output the counts of the QC codes.
        output_qc_codes_counts()

if __name__ == "__main__":
    main()

