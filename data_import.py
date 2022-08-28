import csv
import datetime
from pathlib import Path

from formulas import calc_rh
import archive

def build_row_converter(cols):

    stn_idx, id_idx, date_idx = cols.index('STATION'), cols.index('CALL_SIGN'), cols.index('DATE')
    tmp_idx, dew_idx = cols.index('TMP'), cols.index('DEW')

    def converter(row):
        try:
            stn_num = row[stn_idx]
            stn_id = row[id_idx]
            vt = row[date_idx]
            tmp = row[tmp_idx]
            dew = row[dew_idx]

            if '9999' in tmp or '9999' in dew:
                return None

            stn_num = int(stn_num)

            stn_id = stn_id.strip().lower()
            if len(stn_id) < 4:
                stn_id = 'k' + stn_id

            vt = datetime.datetime.strptime(vt, "%Y-%m-%dT%H:%M:%S")

            tmp = tmp.replace(',', '.')
            dew = dew.replace(',', '.')

            tmp = float(tmp) / 10
            dew = float(dew) / 10

            #tmp = c_to_f(tmp)
            #dew = c_to_f(dew)

            return (stn_num, stn_id, vt, tmp, dew)

        except ValueError:
            return None

    return converter

def round_to_hour(vt):
    hour = vt.hour
    minute = vt.minute
    year, month, day, hour, *_ = vt.timetuple()

    if minute > 30:
        rounded = vt + datetime.timedelta(60 - minute)
        change = 60 - minute
    else:
        rounded = vt - datetime.timedelta(minute)
        change = minute

    year, month, day, hour, *_ = rounded.timetuple()
    rounded = datetime.datetime(year, month, day, hour)

    return (rounded, change)

def aggregate(it):
    current = next(it)
    curr_stn_num, curr_stn_id, curr_vt, curr_tmp, curr_dew = current
    curr_rounded, curr_change = round_to_hour(curr_vt)

    for next_ in it:
        next_stn_num, next_stn_id, next_vt, next_tmp, next_dew = next_
        next_rounded, next_change = round_to_hour(next_vt)

        if next_rounded > curr_rounded:
            yield current

            current = next_
            curr_rounded = next_rounded
            curr_change = next_change
        elif next_change < curr_change:
            current = next_
            curr_rounded = next_rounded
            curr_change = next_change

csv_dir = Path('csv')
csv_files = (x for x in csv_dir.iterdir() if x.is_file() and x.suffix == '.csv')

error_list = []
for fname in csv_files:
    print(f"Loading {fname}")

    with open(fname, 'r') as f, archive.Archive() as ar:
        csvreader = csv.reader(f, delimiter=',', quotechar='"')

        cols = next(csvreader)
        
        stn_idx, id_idx, date_idx = cols.index('STATION'), cols.index('CALL_SIGN'), cols.index('DATE')
        tmp_idx, dew_idx = cols.index('TMP'), cols.index('DEW')

        convert_row = build_row_converter(cols)

        data = (convert_row(r) for r in csvreader)
        data = (r for r in data if r is not None)
        data = aggregate(data)
        data = ((stn_id, vt, tmp, calc_rh(tmp, dew)) for _, stn_id, vt, tmp, dew in data)
        
        for stn, vt, tmp, rh in data:

            # Detect some known errors.
            if stn == '99999' or stn == 'ks80':
                error_list.append((fname, vt, stn))
                continue

            ar.add_update_observation(stn, vt, tmp, rh)

if len(error_list) > 0:
    print("\n\n     ***** Errors *****")
    error_list = tuple(set(error_list))
    for fname, stn in error_list:
        print(fname, stn)

