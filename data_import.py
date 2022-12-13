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

        72794599999: 'kawo', # ARLINGTON MUNICIPAL
        99999904205: 'kawo',
        72794504205: 'kawo',

        72027299999: 'kbvs', # SKAGIT REGIONAL
        99999994282: 'kbvs',
        72027294282: 'kbvs',

        72788594266: 'kclm', # PORT ANGELES FAIRCHILD INTERNATIONAL AIRPORT, WA US
        72788599999: 'kclm',
        99999994266: 'kclm',

        72025400119: 'kcls', # CHEHALIS CENTRALIA AIRPORT, WA US
        94999900119: 'kcls',

        72787399999: 'kcqv', # COLVILLE MUNICIPAL AIRPORT
        72787199999: 'kcqv',

        72781399999: 'kfct', # Vagabond Army Air Field, Yakima

        72798594276: 'kfhr', # FRIDAY HARBOR AIRPORT, WA US
        72798599999: 'kfhr',
        99999994276: 'kfhr',

        74207099999: 'kgrf', # Fort Lewis near Gray AAF

        72784094187: 'khms', # Hanford
        72784099999: 'khms',
        99999994187: 'khms',

        99999924244: 'knej', # SEATTLE NAS

        74201024228: 'knou', # 'know', 'knkw' PORT ANGELES WEATHER BUREAU AIRPORT
        99999924228: 'knou',

        72220804224: 'kors', # EASTSOUND ORCAS ISLAND AIRPORT
        99999904224: 'kors',
        72220899999: 'kors',

        72038800469: 'kplu', # PIERCE CO AIRPORT THUN FIELD

        72788499999: 'krld', # Richland Airport

        99999994248: 'krnt', # RENTON MUNICIPAL

        72793499999: 'krnt', # RENTON MUNICIPAL
        72793494248: 'krnt',

        74201699999: 'ks19', # FRIDAY HARBOR

        72796599999: 'ks88', # ARLINGTON AWOS
        74918599999: 'ks88',

        72792594227: 'kshn', # SHELTON AIRPORT
        72792599999: 'kshn',

        72792699999: 'ktdo', # ED CARLSON MEM FIELD
        72792624241: 'ktdo',
        99999924241: 'ktdo',

        72793894274: 'ktiw', # TACOMA NARROWS AIRPORT, WA US
        72793899999: 'ktiw',
        99999994274: 'ktiw',

        99430099999: 'ktti', # TATOOSH ISLAND WA
        99999924240: 'ktti',

        99999924252: 'kuil', # QUILLAYUTE 

        72791894298: 'kvuo', # VANCOUVER PEARSON AIRPORT, WA US
        72791899999: 'kvuo',
        99999994298: 'kvuo',

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

        99999924134: 'kbno', # Federal Building, Burns, OR

        72698799999: 'kczk', # CASCADE LOCKS STATE
        99999994204: 'kczk',
        72698794204: 'kczk',

        72688399999: 'khri', # Hermiston 
        72688304113: 'khri',
        99999904113: 'khri',

        99999924249: 'klmt', # Klamath

        72688199999: 'kmmv', # MCMINNVILLE MUNICIPAL
        99999994273: 'kmmv',
        72688194273: 'kmmv',

        72695199999: 'konp', # Newport
        99428099999: 'konp',

        72695024285: 'konp', # NEWPORT MUNICIPAL AIRPORT
        72695899999: 'konp',
        72695499999: 'konp',
        72695824285: 'konp',
        72695099999: 'konp',
        99999924285: 'konp',

        72687499999: 'kp88',  # ROME AUTOMATIC METEOROLOGICAL OB

        72695799999: 'kpfc', # PACIFIC CITY STATE

        72687599999: 'kreo', # ROME STATE
        72687594107: 'kreo',
        99999994107: 'kreo',

        74918499999: 'ks33', # CITY CO AIRPORT

        72683604201: 'kspb', # SCAPPOOSE INDUSTRIAL AIRPORT
        99999904201: 'kspb',
        72683699999: 'kspb',

        72020200118: 'ktmk', # TILLAMOOK AIRPORT
        72696399999: 'ktmk',
        94999900118: 'ktmk',
        72020299999: 'ktmk',
        99999924254: 'ktmk',

        72695994281: 'kuao', # AURORA STATE AIRPORT
        72695999999: 'kuao',
        99999994281: 'kuao',

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

        72681604110: 'kjer', # JEROME CO AIRPORT
        99999904110: 'kjer',
        72681699999: 'kjer',

        72073400264: 'kman', # NAMPA MUNICIPAL AIRPORT
        94999900264: 'kman',

        72578624151: 'kmld', # MALAD CITY
        72578699999: 'kmld',
        99999924151: 'kmld',

        72783699999: 'kmlp', # Mullan Pass, this location is down in town.

        99999994182: 'kmyl', # McCall
        72578899999: 'kmyl',
        72578894182: 'kmyl',
        72586499999: 'kmyl',
        72586494182: 'kmyl',

        72587199999: 'kowy', # Middle of nowhere. {42.583} {-116.167} {', US'}

        72578499999: 'kp69', # Powell
        99999904109: 'kp69',
        72578404109: 'kp69',

        72681894194: 'krxe', # REXBURG MADISON CO AIRPORT

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

        72586694178: 'ktwf', # TWIN FALLS SUN VALLEY REGIONAL AIRPORT
        99999994178: 'ktwf',
        72586699999: 'ktwf', 

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

        72677794055: 'kbhk', # Baker
        72677799999: 'kbhk',
        99999994055: 'kbhk',

        72679099999: 'kbtm2', # Butte, but 2.5km away.

        72099899999: 'kem1', # Echo Winifred
        
        72772599999: 'ketm', # Elliston

        99999994010: 'kgsg', #  ST. MARIE

        72768599999: 'kgsg2', #  GLASGOW INDUSTRIAL

        72776099999: 'kgtf', # Great Falls

        72677499999: 'khmm', #  HAMILTON RAVALLI CO
        72075999999: 'khmm',

        99999924035: 'khvr', # Havre Weather Office

        99999994051: 'kjdn', # Jordan
        72768494051: 'kjdn',
        72768499999: 'kjdn',
        
        72100699999: 'klm1', #  LIMA JUDITH GAP

        72679824150: 'klvm', # Livingston
        99999924150: 'klvm',

        72100799999: 'kmm1', #  MIKE MOORE

        72676699999: 'kmqm', # Monida

        72768694017: 'kolf', # Wolf Point
        99999994017: 'kolf',
        72768699999: 'kolf',

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

        72665494054: 'kbyg', # BUFFALO JOHNSON CO AIRPORT
        99999994054: 'kbyg',
        72665499999: 'kbyg',
        72670499999: 'kbyg',

        99999924045: 'kcod2', # Cody but about 12km away, big range in coords for this location.

        72568699999: 'kdgw', # CONVERSE CO
        72568694057: 'kdgw',
        99999994057: 'kdgw',

        72052200443: 'kdub', # DUBOIS MUNICIPAL AIRPORT

        99999900446: 'kecs', # MONDELL FIELD AIRPORT

        72097199999: 'kfwz', # SOUTH PASS
        72097100342: 'kfwz',
        94999900342: 'kfwz',

        72666799999: 'kgey', # South Big Horn Co
        99999924048: 'kgey',
        72666724048: 'kgey',

        72096900341: 'kjpd', # TEN SLEEP
        72096999999: 'kjpd',
        94999900341: 'kjpd',

        72666494173: 'kp60', # Yellowstone Lake

        72034599999: 'kpna', # RALPH WENZ FIELD
        72034594086: 'kpna',

        72576524061: 'kriw', # RIVERTON REGIONAL AIRPORT
        72672099999: 'kriw',
        72576599999: 'kriw',
        72672024061: 'kriw',

        72564899999: 'ksib', # SIBLEY PEAK

        72103600353: 'ktbx', # BOYSEN THERMOPOL
        72103699999: 'ktbx',
        94999900353: 'ktbx',

        72576394053: 'ktor', # TORRINGTON MUNICIPAL AIRPORT
        99999994053: 'ktor',
        72576399999: 'ktor',

        72564699999: 'kvdw', # VEDAUWOO, WY US

        99999900484: 'kw43', # HULETT MUNICIPAL AIRPORT

        72676099999: 'kwey', # West Yellowstone, but different, big range in coords

        72666524062: 'kwrl', # WORLAND MUNICIPAL AIRPORT
        72666599999: 'kwrl',

        72564799999: 'kwtr', # WHITAKER, WY US

        72574399999: 'ocs', # ROCK SPRINGS VORTAC, WY US

        72051699999: 'xafm', # AFTON MUNICIPAL

        99999924016: 'xcas', # CASPER

        72068499999: 'xphi', # PHIFER AFLD

        72051999999: 'xpwl', # POWELL MUNICIPAL

        72051899999: 'kbfr', # Near Fort Bridger

        72467199999: 'lxv', # Near Leadville, CO

        72097799999: 'xxxx', # ALPHA BURNS, WY US

        72074799999: 'xxxx', # ELK MOUNTAIN, WY US

        72061999999: 'xxxx', # FRANCIS E WARREN AFB HELIPORT, WY US

        # North Dakota
        72101899999: 'k1hn', # Hotel Parshall

        72101999999: 'k1in', # India Palermo

        72102299999: 'k1kn', # KILO DONNYBROOK 2

        72102399999: 'k1ln', # LIMA BOWBELLS

        72051800449: 'kbfr', # Fort Bridger

        72529699999: 'kbpp', # Bowman

        72086100287: 'kd50', # Crosby

        72758494038: 'khei', # Hettinger
        72758499999: 'khei',
        99999994038: 'khei',

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

        72651499999: 'kcut', # CUSTER CO
        72651494032: 'kcut',

        99999994056: 'kd07', # Faith Municipal
        72653999999: 'kd07', 
        72653994056: 'kd07',

        72086300289: 'kd60', # Tioga
        94999900289: 'kd60',

        72660699999: 'kefc', # Belle Fourche Municiple

        72651794039: 'kien', # PINE RIDGE AIRPORT
        72651799999: 'kien',
        99999994039: 'kien',

        72683724162: 'kono', # ONTARIO
        72683799999: 'kono',
        99999924162: 'kono',

        72661099999: 'krej', # REDIG SD.

        72660500386: 'kspf', # BLACK HILLS AIRPORT CLYDE ICE FIELD
        72660599999: 'kspf',

        72669199999: 'ky22', # Lemmon
        72669599999: 'ky22',

        # Nebraska
        72563599999: 'kaia', # ALLIANCE MUNICIPAL
        72563524044: 'kaia',
        99999924044: 'kaia',

        72563799999: 'kgrn', # GORDON MUNICIPAL

        72098699999: 'xisn', # INDIA SIDNEY, NE US

        72098299999: 'xepf', # ECHO PINE BLUFFS, NE US

        72098199999: 'xdki', # DELTA KIMBALL, NE US

        72098499999: 'xgsy', # GOLF SIDNEY, NE US

        72098899999: 'xkdx', # KILO DIX, NE US

        72098599999: 'xhgn', # HOTEL GURLEY, NE US
        
        72098399999: 'xfxd', # FOXTROT DIX, NE US

        72097899999: 'xxxx', # BRAVO PINE BLUFFS, NE US

        72097999999: 'xxxx', # CHARLIE GERING, NE US

        # California
        99231099999: 'db23', # Data buoy in the Pacific ocean?
        99238099999: 'db23', 

        99249099999: 'db24', # Data buoy in the Pacific ocean?
        99241099999: 'db24',
        99248099999: 'db24',

        99259099999: 'db25', # Data buoy in the Pacific ocean?
        99258099999: 'db25', 
        99250099999: 'db25', 

        99265099999: 'db226', # Data buoy in the Pacific ocean?
        99260099999: 'db226',
        99264099999: 'db226',

        99294099999: 'db29', # Data buoy in the Pacific ocean?
        99295099999: 'db29',
        99296099999: 'db29',

        72390099999: 'k87q',
        72390093226: 'k87q',

        72480023157: 'kbih',

        72497199999: 'kcic',
        72497393203: 'kcic',
        72497399999: 'kcic',

        72040699999: 'kdvo', # GNOSS FIELD AIRPORT
        72040600135: 'kdvo',

        72389799999: 'kfch', # FRESNO CHANDLER EXECUTIVE
        72389723167: 'kfch', 

        72492723285: 'klvk', # Livermore Municipal AP
        72492723285: 'klvk',
        72492799999: 'klvk',

        72389403181: 'kmmh', # Mammoth Host Springs
        72389499999: 'kmmh', 

        99999903181: 'kmmh', # Near Mammoth Hot Springs

        72290199999: 'kmyf',
        72290303131: 'kmyf',
        72290399999: 'kmyf',
        99999903131: 'kmyf',

        72292603154: 'knfg', 

        72290993115: 'knrs', # Imperial Beach

        72391093111: 'kntd', 

        72286699999: 'ksbd', # SAN BERNARDINO INTERNATIONAL, CA US
        99999923122: 'ksbd', 
        72286623122: 'ksbd', 

        72289793206: 'ksbp', # SAN LUIS OBISPO MCCHESNEY FIELD, CA US
        72289799999: 'ksbp',
        99999993206: 'ksbp',

        72290403178: 'ksdm', # BROWN FIELD MUNICIPAL, CA US
        72290423196: 'ksdm',
        72290499999: 'ksdm',
        99999903178: 'ksdm',

        72290753143: 'ksee', # GILLESPIE FIELD AIRPORT, CA US
        72290799999: 'ksee',
        99999953143: 'ksee',

        72595524259: 'ksiy', # SISKIYOU CO, CA US
        72595599999: 'ksiy',

        72494523293: 'ksjc', # SAN JOSE INTERNATIONAL AIRPORT, CA US

        72297553141: 'ksli', #LOS ALAMITOS ARMY AIR FIELD, CA US 
        72297593106: 'ksli',
        72297599999: 'ksli',
        99999953141: 'ksli',

        72288593197: 'ksmo', # SANTA MONICA MUNICIPAL, CA US
        72288599999: 'ksmo',
        99999993197: 'ksmo',

        72297793184: 'ksna', # SANTA ANA JOHN WAYNE AIRPORT, CA US
        72297799999: 'ksna',

        72394023273: 'ksmx',

        72493893231: 'ksql', # SAN CARLOS, CA US
        72493899999: 'ksql',
        99999993231: 'ksql',

        72495723213: 'ksts', # SANTA ROSA SONOMA CO AIRPORT, CA US
        72495799999: 'ksts',

        72295503174: 'ktoa', # TORRANCE, CA US
        72295503122: 'ktoa',
        72295599999: 'ktoa',
        99999903174: 'ktoa',

        72584793230: 'ktvl', # SOUTH LAKE TAHOE AIRPORT, CA US
        72584799999: 'ktvl',

        74917100479: 'ktsp', # TEHACHAPI MUNICIPAL AIRPORT, CA US

        72584693201: 'ktrk', # TRUCKEE TAHOE, CA US
        72584699999: 'ktrk', # Missing a LOT of data.
        99999993201: 'ktrk',

        72590523275: 'kuki', # UKIAH MUNICIPAL, CA US
        72590599999: 'kuki',

        72482893241: 'kvcb', # NUT TREE, CA US
        99999993241: 'kvcb',
        72482899999: 'kvcb',

        72393099999: 'kvbg', # Near Vandenberg

        72393093214: 'kvbg', # VANDENBERG AFB, CA US
        99999993214: 'kvbg',

        72382523131: 'kvcv', # VICTORVILLE SOUTHERN CALIFORNIA INTERNATIONAL AIRPORT, CA US
        72382599999: 'kvcv',

        72288623130: 'kvny', # VAN NUYS, CA US
        72288699999: 'kvny',

        74505753130: 'kwhp', # LOS ANGELES WHITEMAN AIRPORT, CA US
        74505799999: 'kwhp',

        74505823277: 'kwvi', # WATSONVILLE MUNICIPAL, CA US
        99999923277: 'kwvi',
        74505899999: 'kwvi',

        72393599999: 'kxvw', # VANDENBERG RANGE, CA US

        99729299999: 'ljpc1', # LA JOLLA, CA US

        99817399999: 'mlsc1', # MOSS LANDING S HARBOR, CA US

        74506023239: 'ngz', # ALAMEDA NAS, CA US
        72492523240: 'nrc', # CROWS LANDING, CA US
        69016093114: 'ntk', # TUSTIN MCAF, CA US
        69014093101: 'nzj', # EL TORO MCAS, CA US

        72583799999: 'ksve', # Susanville

        72491699999: 'oar', # MARINA MUNICIPAL
        69007093217: 'oar',
        72491693217: 'oar',


        72286599999: 'ont', # Near Ontario
        72286593180: 'ont', 

        99412099999: 'ptac', # POINT ARENA CA, CA US

        99421099999: 'ptgc', # POINT ARGUELLO CA, US

        72032999999: 'xcab', # CABLE

        72481023203: 'xcaf', # CASTLE AFB, CA US

        72068199999: 'xclv', # CALAVERAS CO MAURY RASMUSSEN FIELD

        72068799999: 'xhem', # HEMET RYAN

        99681099999: 'xhep', # HARVEST EXPERIMENT PLATFORM, US

        99637099999: 'xred', # Redondo Beach, Buoy?

        99640099999: 'xsbs', # SANTA BARB E 12NM SOUTHWEST OF SAN BAR

        74500099999: 'xshw', # SHERIDAN 

        99847599999: 'xmac', # MARTINEZ AMORCO

        99801399999: 'xxxx', # TIJUANA RIVER RESERVE

        99801199999: 'xxxx', # SAN FRANCISCO BAY RESERVE, CA US

        99799799999: 'xxxx', # ELKHORN SLOUGH RESERVE, CA US

        99817599999: 'xxxx', # CAL POLY PIER, CA US

        72095299999: 'xxxx', # GRAY BUTTE FIELD, CA US

        74504099999: 'xxxx', # PILAR POINT AFS CA., US

        99773499999: 'tibc1', 

        # Arizona
        72378199999: 'kgcn', # Grand Canyon
        72378303195: 'kgcn',
        72378399999: 'kgcn',
        99999903195: 'kgcn',

        72278903192: 'ksdl', # SCOTTSDALE, AZ US
        72278999999: 'ksdl',
        99999903192: 'ksdl',

        72375699999: 'ksez', # Sedona
        72375600375: 'ksez',

        72375499999: 'ksjn', # Saint Johns
        72375493027: 'ksjn',
        99999993027: 'ksjn',

        72280023195: 'kyum', # YUMA INTERNATIONAL AIRPORT, AZ US
        72280099999: 'kyum',
        72280199999: 'kyum',

        72378799999: 'lhu', # LAKE HAVASU AWOS, AZ US

        72376099999: 'xbel', # Near FLAGSTAFF

        74916799999: 'xpin', # PINAL AIRPARK

        74917399999: 'xxxx', # COOLIDGE MUNICIPAL

        # Texas
        72270499999: 'kbif', # Fort Biggs Army Air Field

        72265023023: 'maf', # MIDLAND INTERNATIONAL AIRPORT, TX US

        74730099999: 'p07', # SANDERSON TX., US

        # Utah
        72472193025: 'k4bl',
        72472393025: 'k4bl',

        72473599999: 'k4hv', # Near Hanksville AP

        72473523170: 'k4hv', # Hanksville

        72475099999: 'kmlf', # MILFORD MUNICIPAL, UT US
        72479723176: 'kmlf', 
        72475023176: 'kmlf', 

        99999993141: 'kpuc', # Price Carbon Co

        72470093141: 'kpuc', # Price Carbon County
        72470099999: 'kpuc',

        72472599999: 'kt62', # TOOELE, UT US

        72575399999: 'ku16', # EAGLE RANGE, UT US
        99999900480: 'ku16',
        69430099999: 'ku16',
        69133499999: 'ku16',

        99999900481: 'ku19', # GRANITE PEAK FILLMORE AIRPORT, UT US

        72479599999: 'ku24', # DELTA, UT US

        72477199999: 'ku28', # Green River
        72477393132: 'ku28',
        72477399999: 'ku28',
        99999993132: 'ku28',

        74420099999: 'ku67', # ROOSEVELT UT, US

        99999994030: 'kvel', # Near Vernal

        72570594030: 'kvel', # VERNAL MUNICIPAL AIRPORT, UT US
        72570599999: 'kvel',

        72472399999: 'k4bl', # Blanding UT
        72472099999: 'k4bl',

        72473699999: 'u17', # BULLFROG MARINA, UT US

        72056999999: 'xknb', # KANAB MUNICIPAL

        72057299999: 'xbof', # BOLINDER FIELD TOOELE VALLEY

        # Nevada
        72582794190: 'kawh', # Wild horse reservoir
        72582799999: 'kawh',
        99999994190: 'kawh',

        74614099999: 'kins', # Near Indian Springs

        74614123141: 'kins', # Indian Springs
        74614199999: 'kins',

        72386599999: 'klsv', # Near Nellis AFB

        72386523112: 'klsv', # Nellis AFB
        74614023112: 'klsv',

        72484653123: 'kvgt', # LAS VEGAS AIR TERMINAL, NV US
        72484699999: 'kvgt',
        99999953123: 'kvgt',

        72484499999: 'ktnx', # TONOPAH TEST RANGE, NV US

        72582699999: 'ku31', # AUSTIN, NV US

        72582499999: 'p68', # EUREKA RAMOS, NV US

        72055199999: 'xmta', # MINDEN TAHOE AIRPORT, NV US

        72012999999: 'xaan', # AUSTIN AIRPORT, NV US

        72582099999: 'xelk', # ELKO, NV US

        72485699999: 'xxxx', # HAWTHORNE INDUSTRIAL, NV US

        # New Mexico
        69001499999: '2c2', # Near White Sands

        72362099999: 'konm',
        72362093040: 'konm',

        72365623049: 'ksaf', # SANTA FE CO MUNICIPAL, NM US
        72365699999: 'ksaf',

        72272093063: 'ksvc', # Grant CO. Silver City
        72272193063: 'ksvc',

        72367623048: 'ktcc', # TUCUMCARI MUNICIPAL, NM US
        72367699999: 'ktcc',

        72365399999: 'kzab', # ALBUQUERQUE RADAR, NM US

        74637599999: 'p70', # CLINES CORNER, NM US

        72364099999: 'xstj', # DONA ANA CO INTERNATIONAL JETPORT ARPT SANTA TERESA

        72096399999: 'xxxx', # MCGREGOR RANGE BASE CAMP, NM US

        74638099999: 'xxxx', # MELROSE GUNNERY RANGE, NM US

        # Colorado
        72053699999: 'grby', # GRANBY GRAND CO, CO US

        72469023062: 'kden', # Near DIA, 

        72053499999: 'keik', # ERIE MUNICIPAL, CO US
        72053400161: 'keik', 

        99999993007: 'kguc', # Near Gunnison

        72467894033: 'ksbs', # STEAMBOAT SPRINGS ADAMS FIELD, CO US
        72467899999: 'ksbs',

        72464603028: 'kspd', # SPRINGFIELD COMANCHE NATIONAL GRASSLAND, CO US
        72464699999: 'kspd',
        99999903028: 'kspd',

        72054400168: 'kstk', # STERLING MUNICIPAL AIRPORT, CO US

        74000203042: 'kvtp', # LA VETA PASS, CO US
        74000299999: 'kvtp',
        99999903042: 'kvtp',

        72468499999: 's29', # SALIDA, CO US

        72639100385: 'xbmc', # BALD MOUNTAIN COTTONWOOD PASS, CO US

        72429300376: 'xmtw', # MOUNT WERNER, CO US

        72639200424: 'xsun', # SUNLIGHT, CO US

        72639600422: 'xwkp', #  WILKERSON PASS, CO US

        72098999999: 'xxxx', # LIMA GROVER, CO US

        72099299999: 'xxxx', # NOVEMBER GROVER, CO US

        99999994044: 'xxxx', #

        72098799999: 'xxxx', # JULIET PEETZ, CO US

        72084499999: 'xxxx', # SPANISH PEAKS, CO US

        72042499999: 'xxxx', # CAMP RED DEVIL, CO US

        72099399999: 'xxxx', # OSCAR GROVER, CO US

        72099199999: 'xxxx', # MIKE HAXTUN, CO US

        72463099999: 'xxxx', # LAMAR, CO US

        72064199999: 'xxxx', # SCHRIEVER AFB, CO US

        72064299999: 'xxxx', # CHEYENNE MOUNTAIN AS, CO US

        # Mexico
        76075399999: 'cjs', # ABRAHAM GONZALEZ INTERNATIONAL, MX
        76075199999: 'cjs',

        76113099999: 'masm', # ALTAR SON.

        76122099999: 'mcgs', # Casas Grandes Chihuahua

        76050099999: 'mebc', # ENSENADA BC, MX

        76005399999: 'mmml', # GENERAL RODOLFO SANCHEZ TABOADA INTERNATIONAL
        76005499999: 'mmml',
        76005199999: 'mmml',

        76001499999: 'mmtj', # TIJUANA, MX

        76118099999: 'mndg', # PILARES DE NACOZARI SON., MX

        76040099999: 'mnlb', # Baha, Nuevo Leon

        76055099999: 'msfb', # SAN FELIPE BCN, MX

        76061099999: 'mpps', # PUERTO PENASCO SON., MX

        76001199999: 'tij', # TIJUANA

        76062099999: 'xppm', # PUERTO PENASCO, MX

        76001399999: 'xxxx', # GENERAL ABELARDO L RODRIGUEZ INTERNATIONAL, MX

        # Lots of sites I'm skipping because of too few values and I didn't want to organize them
        # into states.
        72595999999: '1o5',
        72571699999: '1v1',
        69031003125: '1y7',
        72468699999: '2v9',
        74619099999: '4su',
        72462699999: '4v5',
        72476699999: '6v8',
        72467999999: 'c96',
        99999993134: 'cqt',
        99229099999: 'db22',
        99274099999: 'db27',
        99281099999: 'db28',
        99647099999: 'db64',
        72365999999: 'e23',
        74733099999: 'e28',
        72274399999: 'emx',
        72469799999: 'fcl',
        72375799999: 'fsx',

        72233303069: 'k04v',
        72233399999: 'k04v',

        99999900170: 'k05u',
        72054699999: 'k05u',

        72038500419: 'k0co',

        99999903033: 'k0e0',
        72077299999: 'k0e0',

        72374593139: 'k0e4',
        72272699999: 'k13a',
        72052499999: 'k1v6',

        72026294076: 'k20v',
        72026299999: 'k20v',
        99999994076: 'k20v',
        99999994076: 'k20v',

        72269099999: 'k2c2',

        99999900319: 'k2v5',
        72093499999: 'k2v5',

        99999900157: 'k33v',

        99999904133: 'k36u',
        72056599999: 'k36u',

        72004699999: 'k3a6',
        72005999999: 'k40g',

        72268599999: 'k4cr',
        72365599999: 'k4my',

        72365703014: 'k4sl',
        72365799999: 'k4sl',
        99999903014: 'k4sl',

        99999903053: 'k5t6',

        72261603032: 'k6r6',
        72261603032: 'k6r6',
        72261699999: 'k6r6',
        99999903032: 'k6r6',
        99999903032: 'k6r6',

        72484599999: 'k9bb',

        72317153144: 'k9l2',
        72317199999: 'k9l2',

        72595894299: 'kaat',
        72595894299: 'kaat',
        99999994299: 'kaat',
        99999994299: 'kaat',

        72281703068: 'kabh',

        72364703034: 'kaeg',
        99999903034: 'kaeg',
        99999903034: 'kaeg',
        72364799999: 'kaeg',

        72052803064: 'kaej',

        74531093065: 'kaff',
        74531099999: 'kaff',

        72052900429: 'kaib',
        72033353175: 'kajo',
        72053100158: 'kajz',

        72281353146: 'kalk',
        72281353146: 'kalk',

        72053200159: 'kank',

        72466693067: 'kapa',
        72466693067: 'kapa',
        72466699999: 'kapa',
        99999993067: 'kapa',
        99999993067: 'kapa',

        72495593227: 'kapc',
        72495593227: 'kapc',
        72495599999: 'kapc',

        72564399999: 'karl',

        72467693073: 'kase',
        72467693073: 'kase',
        72467699999: 'kase',
        99999993073: 'kase',
        99999993073: 'kase',

        72026723224: 'kaun',
        72026799999: 'kaun',
        99999923224: 'kaun',
        99999923224: 'kaun',

        72292023191: 'kavx',
        72292023191: 'kavx',
        72292099999: 'kavx',

        72041100137: 'kaxx',

        72583524119: 'kb23',
        72583599999: 'kb23',

        72483793216: 'kbab',
        72483799999: 'kbab',

        72065200433: 'kban',
        72065299999: 'kban',

        72475623159: 'kbce',
        72475623159: 'kbce',
        72475699999: 'kbce',

        72053300160: 'kbdu',

        72469903065: 'kbjc',
        72469999999: 'kbjc',
        99999903065: 'kbjc',
        99999903065: 'kbjc',

        72282553145: 'kbjn',
        72282599999: 'kbjn',

        72584523225: 'kblu',
        72584523225: 'kblu',

        72056724180: 'kbmc',
        72595699999: 'kbny',

        72564999999: 'kbrx',

        72286723156: 'kbuo',
        72286799999: 'kbuo',

        72074100269: 'kbvu',

        72064400226: 'kbxk',
        72064499999: 'kbxk',

        74611003182: 'kbys',
        74611099999: 'kbys',

        72495023254: 'kccr',
        72493623254: 'kccr',

        72206103038: 'kccu',
        72206199999: 'kccu',
        99999903038: 'kccu',

        72476793069: 'kcez',
        72476799999: 'kcez',
        99999993069: 'kcez',

        72274953128: 'kchd',
        72274999999: 'kchd',
        99999953128: 'kchd',

        72063500221: 'kcmr',

        72289903179: 'kcno',
        72289999999: 'kcno',
        99999903179: 'kcno',

        72477693075: 'kcny',
        72477699999: 'kcny',
        99999993075: 'kcny',

        72210103039: 'kcpw',
        72210199999: 'kcpw',
        99999903039: 'kcpw',

        72267703027: 'kcqc',
        99999903027: 'kcqc',
        72267799999: 'kcqc',

        72287493134: 'kcqt',

        72292703177: 'kcrq',
        72292703177: 'kcrq',
        72292799999: 'kcrq',
        99999903177: 'kcrq',
        99999903177: 'kcrq',

        74917900392: 'kcvh',

        72268623008: 'kcvs',
        72268699999: 'kcvs',

        72054900171: 'kcxp',

        74718603164: 'kczz',
        74718699999: 'kczz',
        99999903164: 'kczz',
        99999903164: 'kczz',

        72469099999: 'kdnr',

        74003024103: 'kdpg',
        74003099999: 'kdpg',
        99999924103: 'kdpg',
        99999924103: 'kdpg',

        72462593005: 'kdro',
        72462593005: 'kdro',
        72462599999: 'kdro',
        99999993005: 'kdro',
        99999993005: 'kdro',

        72221403070: 'kdux',
        72221499999: 'kdux',

        72222100444: 'kdwx',
        72092299999: 'kdxz',

        99999900349: 'ke11',
        72103299999: 'ke11',

        74948400395: 'ke16',

        72367593092: 'ke33',
        72367599999: 'ke33',
        99999993092: 'ke33',

        72015103049: 'ke38',
        99999903049: 'ke38',
        99999903049: 'ke38',
        72015199999: 'ke38',

        72274799999: 'ke74',

        99999900258: 'ke80',
        72072499999: 'ke80',

        72057600174: 'kedu',

        72594024213: 'keka',
        72594099999: 'keka',

        72051724165: 'kemm',

        74704303165: 'kemt',
        74704399999: 'kemt',
        99999903165: 'kemt',

        72577504111: 'kevw',
        72577504111: 'kevw',
        72577599999: 'kevw',

        72468094015: 'kfcs',
        72468099999: 'kfcs',
        99999994015: 'kfcs',
        99999994015: 'kfcs',

        72278303185: 'kffz',
        72278303185: 'kffz',
        72278399999: 'kffz',
        99999903185: 'kffz',
        99999903185: 'kffz',

        72085200280: 'kfly',
        72053500162: 'kfmm',
        74948500396: 'kfot',

        72261823091: 'kfst',
        72261823091: 'kfst',
        72265499999: 'kfst',
        72261899999: 'kfst',
        99999923091: 'kfst',
        99999923091: 'kfst',
        72265423091: 'kfst',

        72469400450: 'kftg',

        72297603166: 'kful',
        72297603166: 'kful',
        72297699999: 'kful',
        99999903166: 'kful',
        99999903166: 'kful',

        74724003148: 'kgbn',
        74724023168: 'kgbn',
        74724099999: 'kgbn',

        72262023055: 'kgdp',
        72262023055: 'kgdp',

        72278753126: 'kgeu',
        72278799999: 'kgeu',
        99999953126: 'kgeu',
        99999953126: 'kgeu',

        72210703056: 'kgnc',
        99999903056: 'kgnc',
        99999903056: 'kgnc',
        72210799999: 'kgnc',

        72362593057: 'kgnt',
        72362599999: 'kgnt',

        74948600397: 'kgoo',

        72278803186: 'kgyr',
        72278899999: 'kgyr',
        99999903186: 'kgyr',

        72064600228: 'khaf',

        72038799999: 'khbb',

        72053700163: 'kheq',

        69002093218: 'khgt',
        69002099999: 'khgt',

        72295603167: 'khhr',
        72295603167: 'khhr',
        72295699999: 'khhr',
        99999903167: 'khhr',
        99999903167: 'khhr',

        72389853119: 'khjo',
        72389853119: 'khjo',
        99999953119: 'khjo',
        99999953119: 'khjo',
        72389899999: 'khjo',

        74732023002: 'khmn',
        74732099999: 'khmn',
        99999923002: 'khmn',
        99999923002: 'khmn',

        72209653127: 'khnd',
        72209699999: 'khnd',
        99999953127: 'khnd',
        99999953127: 'khnd',

        72268893034: 'khob',
        72268899999: 'khob',
        99999993034: 'khob',
        99999993034: 'khob',

        72039300129: 'khrx',

        72473323170: 'khve',

        72585093228: 'khwd',
        72585093228: 'khwd',
        72493593228: 'khwd',
        72493599999: 'khwd',

        72566594073: 'kibm',
        72566599999: 'kibm',
        99999994073: 'kibm',
        99999994073: 'kibm',

        72378853135: 'kifp',
        72378899999: 'kifp',
        99999953135: 'kifp',
        99999953135: 'kifp',

        72265623040: 'kink',
        72265623040: 'kink',

        72468903026: 'kitr',
        72468903026: 'kitr',
        99999903026: 'kitr',
        99999903026: 'kitr',
        72468999999: 'kitr',

        72278623104: 'kiwa',
        72278699999: 'kiwa',

        72382693194: 'kiyk',
        72382699999: 'kiyk',

        72376293244: 'kiza',
        99999993244: 'kiza',
        72376299999: 'kiza',

        72016599999: 'kl35',

        72067599999: 'kl63',
        69017099999: 'kl63',

        72463603013: 'klaa',
        72463603013: 'klaa',
        72463699999: 'klaa',
        99999903013: 'klaa',
        99999903013: 'klaa',

        72479694128: 'klgu',
        72479694128: 'klgu',
        72479699999: 'klgu',

        72061400205: 'klhm',


        72053800164: 'klmo',

        72097600345: 'klsb',

        72278523111: 'kluf',
        72278599999: 'kluf',

        72367723054: 'klvs',
        72367723054: 'klvs',

        72467393009: 'klxv',
        72467393009: 'klxv',
        72467399999: 'klxv',
        99999993009: 'klxv',
        99999993009: 'klxv',

        74504693242: 'kmae',
        74504693242: 'kmae',
        99999993242: 'kmae',
        74504699999: 'kmae',

        72201899999: 'kmai',

        72483623208: 'kmcc',
        72483699999: 'kmcc',

        72481523257: 'kmce',
        72481523257: 'kmce',
        72481599999: 'kmce',
        99999923257: 'kmce',
        99999923257: 'kmce',

        72232403071: 'kmdd',
        72232499999: 'kmdd',

        72481099999: 'kmer',

        72295303183: 'kmhv',
        72295399999: 'kmhv',

        72215503040: 'kmnh',
        99999903040: 'kmnh',
        99999903040: 'kmnh',
        72215599999: 'kmnh',

        72492623258: 'kmod',
        72492623258: 'kmod',

        72491523259: 'kmry',
        72491523259: 'kmry',
        72491523245: 'kmry',
        99999923259: 'kmry',
        72491599999: 'kmry',

        72476593013: 'kmtj',
        72476593013: 'kmtj',
        72476599999: 'kmtj',
        99999993013: 'kmtj',
        99999993013: 'kmtj',

        72289093136: 'kmws',
        72289099999: 'kmws',

        72220303041: 'kmyp',
        99999903041: 'kmyp',
        99999903041: 'kmyp',
        72220399999: 'kmyp',

        72483893205: 'kmyv',
        72483893205: 'kmyv',
        72483899999: 'kmyv',
        99999993205: 'kmyv',
        99999993205: 'kmyv',

        74612093104: 'knid',
        74612093104: 'knid',

        72281023199: 'knjk',
        72281023199: 'knjk',

        74702023110: 'knlc',

        74505499999: 'knnz',

        72291093116: 'knsi',
        72291093116: 'knsi',
        72291099999: 'knsi',

        72292593117: 'knuc',
        72292599999: 'knuc',

        74509023244: 'knuq',
        74509099999: 'knuq',

        72292800369: 'knxf',

        69015093121: 'knxp',
        69015093121: 'knxp',

        72061500206: 'ko22',
        72019300117: 'ko54',
        72590099999: 'ko64',
        72093599999: 'ko69',
        72102799999: 'ko86',
        72584899999: 'ko87',

        72264803031: 'kodo',
        72264803031: 'kodo',
        72264899999: 'kodo',
        99999903031: 'kodo',
        99999903031: 'kodo',

        72575024126: 'kogd',
        72575024126: 'kogd',
        72575099999: 'kogd',

        72293453121: 'kokb',
        72293453121: 'kokb',
        72293499999: 'kokb',
        99999953121: 'kokb',
        99999953121: 'kokb',

        72272803196: 'kols',
        72272803196: 'kols',
        72272899999: 'kols',
        99999903196: 'kols',
        99999903196: 'kols',

        74704003102: 'kont',
        74704003102: 'kont',
        74704099999: 'kont',

        74504893210: 'kove',
        74504893210: 'kove',
        99999993210: 'kove',
        99999993210: 'kove',
        74504899999: 'kove',

        72392793110: 'koxr',
        72392793110: 'koxr',
        72392799999: 'koxr',

        72487099999: 'kp38',
        72487003163: 'kp38',

        72477003170: 'kp68',
        72582403170: 'kp68',
        72477099999: 'kp68',
        99999903170: 'kp68',
        99999903170: 'kp68',

        72374500374: 'kpan',

        72493723289: 'kpao',
        72493799999: 'kpao',
        99999923289: 'kpao',

        72216903058: 'kpeq',
        72216999999: 'kpeq',
        99999903058: 'kpeq',

        72382023182: 'kpmd',
        72382023182: 'kpmd',
        72382099999: 'kpmd',

        72288703180: 'kpoc',
        72288799999: 'kpoc',
        99999903180: 'kpoc',

        72053900165: 'kpso',

        72286893138: 'kpsp',
        72286893138: 'kpsp',
        72286899999: 'kpsp',
        99999993138: 'kpsp',
        99999993138: 'kpsp',

        72564499999: 'kpum',

        72064500227: 'kpvf',

        69236499999: 'kqav',

        99999904134: 'kqca',
        69258499999: 'kqcr',
        69058499999: 'kqcx',
        69141499999: 'kqdb',
        69811499999: 'kqgy',
        69307499999: 'kqhg',
        69756499999: 'kqiz',
        69249499999: 'kqmp',
        69269499999: 'kqmy',

        69541499999: 'kqnl',

        69542499999: 'kqnm',

        69633499999: 'kqog',

        69645499999: 'kqos',

        69646499999: 'kqot',

        69841499999: 'kqxg',

        72286903171: 'kral',
        72286999999: 'kral',
        99999903171: 'kral',

        72267523021: 'kree',

        72494693232: 'krhv',
        72494699999: 'krhv',
        99999993232: 'krhv',

        72571703016: 'kril',
        72571799999: 'kril',
        99999903016: 'kril',

        74505653120: 'krnm',
        74505699999: 'krnm',
        99999953120: 'krnm',

        72276403029: 'krqe',
        72276499999: 'krqe',
        99999903029: 'krqe',

        72267823052: 'krtn',
        72267899999: 'krtn',

        72083900279: 'krts',

        72669024057: 'krwl',
        72574524057: 'krwl',
        72669099999: 'krwl',

        72033900121: 'kryn',

        72595899999: 'ks11',

        72052100475: 'ksaa',

        72274793084: 'ksad',
        99999993084: 'ksad',

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

        72793599999: 'kbfi', # Seattle Boeing Field
        72793524234: 'kbfi',
        99999924234: 'kbfi',

        72797624217: 'kbli', # BELLINGHAM INTERNATIONAL AIRPORT

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

        72797524255: 'knuw', # WHIDBEY ISLAND NAS
        69023024255: 'knuw',
        72074924255: 'knuw',
        99999924255: 'knuw',

        72792024227: 'kolm', # OLYMPIA AIRPORT, WA US
        99999924227: 'kolm',

        72789099999: 'komk', # Omak
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

        72793024233: 'ksea', # SEATTLE TACOMA AIRPORT, WA US
        99999924233: 'ksea',

        72785699999: 'ksff', # Felts Field
        99999994176: 'ksff',
        72785694176: 'ksff',
        72785694176: 'ksff',

        72785524114: 'kska', # Fairchild AFB
        72785599999: 'kska',
        69148499999: 'kska',

        72781524237: 'ksmp', # Stampede Pass
        99999924237: 'ksmp',

        74206024207: 'ktcm', # TACOMA MCCHORD AFB, WA US
        74206099999: 'ktcm',
        99999924207: 'ktcm', # TACOMA MCCHORD AFB

        72797094240: 'kuil', # QUILLAYUTE AIRPORT
        99999994240: 'kuil',

        72781024243: 'kykm', # Yakima Airport
        99999924243: 'kykm',

        # Oregon
        72791094224: 'kast', # ASTORIA AIRPORT PORT OF
        99999994224: 'kast',
        99999924246: 'kast',

        72688624130: 'kbke', # Baker City
        99999924130: 'kbke',

        72683094185: 'kbno', # BURNS MUNICIPAL AIRPORT
        72683894185: 'kbno',
        72683099999: 'kbno',

        72598524267: 'kbok', # BROOKINGS STATE AIRPORT
        72598599999: 'kbok',
        72598199999: 'kbok',
        72598099999: 'kbok',
        72036524267: 'kbok',
        72036599999: 'kbok',

        72694599999: 'kcvo', # CORVALLIS MUNICIPAL
        99999924248: 'kcvo',
        72694524202: 'kcvo', # CORVALLIS MUNICIPAL
        99999924202: 'kcvo',

        72693024221: 'keug',  # EUGENE MAHLON SWEET AIRPORT
        99999924221: 'keug',

        72698694261: 'khio', # PORTLAND HILLSBORO AIRPORT, OR US
        72698699999: 'khio',
        99999994261: 'khio',

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

        72694024232: 'ksle', # SALEM AIRPORT MCNARY FIELD, OR US
        99999924232: 'ksle',

        72698599999: 'kttd', # PORTLAND TROUTDALE
        99999924242: 'kttd',
        72698524242: 'kttd',

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

        72578024156: 'kpih', # POCATELLO REGIONAL AIRPORT
        99999924156: 'kpih',

        72686599999: 'ksmn', # Salmon
        72686199999: 'ksmn',
        72686524196: 'ksmn',

        # Montana
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

        72677624036: 'klwt', # Lewistown
        99999924036: 'klwt',

        74230024037: 'kmls', # Miles City
        72667524037: 'kmls',
        99999924037: 'kmls',

        72773024153: 'kmso', # Missoula
        99999924153: 'kmso',

        72768794028: 'ksdy', # Sidney
        72768799999: 'ksdy',
        99999994028: 'ksdy',

        72679124135: 'kbtm', # Butte
        72678524135: 'kbtm',
        72774024135: 'kbtm',
        99999924135: 'kbtm',

        # Wyoming
        72671024164: 'kbpi', # BIG PINEY MARBLETON AIRPORT
        72671099999: 'kbpi',

        72670099999: 'kcod', # Cody
        72670024045: 'kcod',
        72666699999: 'kcod',

        72569024089: 'kcpr', # CASPER NATRONA CO INTERNATIONAL AIRPORT
        99999924089: 'kcpr',

        72564024018: 'kcys', # Cheyenne

        72665099999: 'kgcc', # GILLETTE GILLETTE C
        99999994023: 'kgcc',
        72665094023: 'kgcc',
        72666299999: 'kgcc',
        72665599999: 'kgcc',

        72577624166: 'kjac', # JACKSON HOLE AIRPORT
        99999924166: 'kjac',
        72577699999: 'kjac',

        72576024021: 'klnd', # LANDER HUNT FIELD AIRPORT
        99999924021: 'klnd',

        72574124027: 'krks', # Rock Springs
        72574424027: 'krks',

        72666024029: 'kshr', # Sheridan
        99999924029: 'kshr',

        # North Dakota
        72764524012: 'kdik', # Dickinson
        72763024012: 'kdik',
        99999924012: 'kdik',

        72767094014: 'kisn', # Williston
        99999994014: 'kisn', 

        # South Dakota
        72662794037: 'k2wx', # Buffalo
        72662799999: 'k2wx',

        72662024090: 'krap', # RAPID CITY REGIONAL AIRPORT
        99999924090: 'krap',

        72662599999: 'krca', # ELLSWORTH AFB
        72662524006: 'krca',
        99999924026: 'krca',

        72597524235: 'ksxt', # SEXTON SUMMIT
        99999924235: 'ksxt',

        # Nebraska !
        72566024028: 'kbff',

        72563624017: 'kcdr', # CHADRON MUNICIPAL AIRPORT
        72563699999: 'kcdr',
        99999924017: 'kcdr',

        # Colorado
        72462023061: 'kals', 

        72466093037: 'kcos', # Colorado Springs

        72565003017: 'kden', # Denver IAP
        72467003017: 'kden',

        72467523063: 'kege', # EAGLE CO AIRPORT, CO US
        72467599999: 'kege',

        72476023066: 'kgjt', # Grand Junction

        72467799999: 'kguc', # Gunnison
        72467793007: 'kguc',

        72464093058: 'kpub', 

        72462703011: 'ktex', # TELLURIDE REGIONAL AIRPORT, CO US
        72462799999: 'ktex',
        99999903011: 'ktex',

        # California
        72594524283: 'kacv',

        72384023155: 'kbfl',

        72381523161: 'kdag', # Barstow

        72389093193: 'kfat', # Fresno

        72295023174: 'klax', 

        72297023129: 'klgb',

        72289523243: 'klpc', # LOMPOC AIRPORT, CA US
        72289599999: 'klpc',
        99999923243: 'klpc',

        72483523206: 'kmhr',
        72483323206: 'kmhr',
        72483399999: 'kmhr',

        72293193107: 'knkx', 
        72293199999: 'knkx',
        72293093107: 'knkx',

        72290693112: 'knzy', 

        72389523149: 'kptv', # PORTERVILLE MUNICIPAL AIRPORT, CA US
        72389599999: 'kptv',
        99999923149: 'kptv',

        72592024257: 'krdd', # Redding 
        72591524257: 'krdd',

        72492023237: 'ksck', # Stockton

        72383023187: 'ksdb', # SANDBERG, CA US

        72494023234: 'ksfo', # San Francisco 

        4504099999: 'xpps', # PILAR POINT AFS

        # Arizona
        72274803914: 'kcgz',
        72274899999: 'kcgz',
        99999903914: 'kcgz',

        72272099999: 'kdug', # Douglas Bisbee IAP
        72273593026: 'kdug', 
        72273593026: 'kdug', 

        72273003124: 'kfhu', # FORT HUACHUCA, AZ US
        72273099999: 'kfhu',

        72375503103: 'kflg', # Flagstaff
        72375003103: 'kflg', 
        72375099999: 'kflg', 

        72374023194: 'kinw',
        
        72278023183: 'kphx', # Phoenix

        72372323184: 'kprc', # Prescott
        72372123184: 'kprc',

        72374703101: 'ksow', # SHOW LOW AIRPORT, AZ US
        72374799999: 'ksow',
        99999903101: 'ksow',

        72274023160: 'ktus', # Tucson
        72274023160: 'ktus',

        # New Mexico
        72365023050: 'kabq',

        72269393097: 'kalm', # White Sands
        72269399999: 'kalm',
        99999993097: 'kalm',
        99999993097: 'kalm',

        72267603035: 'kats', # ARTESIA MUNICIPAL AIRPORT, NM US
        72267699999: 'kats',
        99999903035: 'kats',
        99999903035: 'kats',

        72360023051: 'kcao',

        72268923077: 'kcvn', # CLOVIS MUNICIPAL AIRPORT, NM US
        72268999999: 'kcvn',
        99999923077: 'kcvn',
        99999923077: 'kcvn',

        72362723081: 'kgup', 

        72365493091: 'klam', # Los Alamos
        72365499999: 'klam',

        72269593041: 'klru', # LAS CRUCES MUNICIPAL AIRPORT, NM US
        72269599999: 'klru',
        99999993041: 'klru',

        72268023009: 'krow', # Roswell

        72268393083: 'ksrr', # RUIDOSO SIERRA BLANCA AIRPORT, NM US
        72268399999: 'ksrr',
        99999993083: 'ksrr',

        # Texas
        72270023044: 'kelp', # El Paso

        72262699999: 'kmrf', # Marfa AP
        72264093035: 'kmrf',
        72264099999: 'kmrf',

        # Utah
        72475593129: 'kcdc',

        72572424174: 'kpvu', # PROVO AIRPORT, UT US
        72572499999: 'kpvu',
        99999924174: 'kpvu',

        72475423186: 'ksgu', # ST GEORGE MUNICIPAL AIRPORT, UT US
        72475499999: 'ksgu',
        99999923186: 'ksgu',

        72572024127: 'kslc',

        # Nevada
        72486023154: 'kely', 

        72386023169: 'klas',

        72580524172: 'klol', # Lovelock

        72488023185: 'krno', # Reno

        72485523153: 'ktph', # Tonapah
        72485523153: 'ktph', 

        72583024128: 'kwmc', # Winnemuca 

        # Mexico!
}

skip_stations_counts = {}

def output_skip_stations_counts():

    # Organize by site
    items_by_site = {}
    for stn_num, (count, earliest, latest) in skip_stations_counts.items():
        site = skip_stations[stn_num]
        if site not in items_by_site:
            items_by_site[site] = (count, earliest, latest)
        else:
            i_count, i_earliest, i_latest = items_by_site[site]
            count += i_count
            if i_earliest < earliest:
                earliest = i_earliest
            if i_latest > latest:
                latest = i_latest

            items_by_site[site] = (count, earliest, latest)

    items = [(site, count, earliest, latest) for site, (count, earliest, latest) 
            in items_by_site.items()]

    items.sort(reverse=True, key=lambda x: x[1])

    log_file.write("\n\nCounts of skipped stations in the data.\n")
    for item in items:
        site, count, earliest, latest = item

        days = count / 24.0
        years = days / 365.25

        years_covered = (latest - earliest).total_seconds() / 60.0 / 60.0 / 24.0 / 365.25 

        log_file.write(f"{site:5s} {years:.0f} years of data covering {years_covered:.0f} years" \
                + f" from {earliest} to {latest}\n")
    
###################################################################################################
# Stations not configured in stations or skip_stations
###################################################################################################
stations_not_configured = {}

def output_stations_not_configured():
    log_file.write("\n\nStations not configured or skipped.\n")
    for stn_num, (site_id, name) in stations_not_configured.items():
        log_file.write(f"{stn_num} - {site_id} - {name}\n")

###################################################################################################
# Check that there's no overlap in the dictionaries
###################################################################################################
def check_dictionaries():
    duplicates = []

    for key in stations:
        if key in skip_stations:
            duplicates.append(key)

    if len(duplicates) > 0:

        print(f"\n{len(duplicates)} duplicate stations found in both stations and skip_stations.\n")
        log_file.write(f"\n{len(duplicates)} duplicate stations found in both stations and skip_stations.\n")

        for stn in duplicates:
            print(stn)
            log_file.write(f"   {stn}\n")

        return False

    return True

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

            # Keep track of how many observations are in the skip stations.
            if stn_num in skip_stations:
                if stn_num not in skip_stations_counts:
                    skip_stations_counts[stn_num] = (1, vt, vt)
                else:
                    count, earliest, latest = skip_stations_counts[stn_num]
                    count += 1
                    if vt < earliest:
                        earliest = vt
                    if vt > latest:
                        latest = vt
                    skip_stations_counts[stn_num] = (count, earliest, latest)

            # Check if we've manually renamed the site_id for this station number
            if stn_num in stations:
                stn_id_forced = stations[stn_num]
                if stn_id_forced != stn_id:
                    name = f"{stn_id} {name}"
                    stn_id = stn_id_forced

            # Check if this station is not configured and is only automatically handled.
            if stn_num not in stations and stn_num not in skip_stations and \
                    stn_num not in stations_not_configured:
                stations_not_configured[stn_num] = (stn_id, name)

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

    if not check_dictionaries():
        return

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

        # Output number of observationsfor skipped stations.
        output_skip_stations_counts()

        # Output stations not configured in stations or skip_stations lists.
        output_stations_not_configured()

if __name__ == "__main__":
    main()

