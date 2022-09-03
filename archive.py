"""An archive of point observations.

 This file is a library to serve as an interface to an sqlite3 file storing observation data.

 Client programs will download and import observations, or query that data for
 plotting or verification.

"""
from datetime import datetime, date, timedelta
from pathlib import Path
import sqlite3

class Archive:
    """An archive of observations.

    All data is point focused, meaning it is associated with an 
    observation site.
    """
    def __init__(self, root=None):
        """Create or open an archive.

        If the root argument is specified it will try to open or connect
        to a sqlite3 database in that directory. If no archive is present,
        one will be created to include creating the path if it does not
        exist. If root is None it will try to create an archive rooted 
        in the same path as this file.

        Keyword Arguments:
        root -- the directory where the database files are/will be stored.
        """
        DB_FILE = 'obs.db'

        # Check for and make the archive folder
        if root is None:
            root = Path(__file__).parent / "default_archive"
        else:
            root = Path(root)
        root.resolve()
        root.mkdir(parents=True, exist_ok=True)
        print("Archive root is", root)

        # Check if the databases exists.
        db_file = root / DB_FILE
        if not db_file.exists():
            print("No database found, one will be created.")

        self.db_conn = sqlite3.connect(db_file,
                                       detect_types=sqlite3.PARSE_DECLTYPES
                                       | sqlite3.PARSE_COLNAMES)
        print("Connected to", db_file)

        # Set up the tables.
        c = self.db_conn.cursor()
        c.executescript("""
                CREATE TABLE IF NOT EXISTS obs(
                    site_id     TEXT      NOT NULL, -- e.g. kmso, kgpi
                    valid_time  TIMESTAMP NOT NULL,
                    temperature REAL      NOT NULL, -- Temperature in C.
                    dew         REAL      NOT NULL, -- Dew Point in C.
                    pres        REAL,               -- Pressure in hPa
                    UNIQUE (site_id, valid_time)
                );
            """)

        return

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.db_conn.commit()
        self.db_conn.close()

    def add_update_observation(self, site, valid_time, temperature, dew_point, pressure=None):
        """Add an observation to the archive.

        The valid_time of the observation is rounded to the nearest hour
        so that it will match up properly with forecasts.

        Arguments:
        site -- a string with the site id: e.g. kmso. Case insensitive.
        valid_time -- a datetime.datetime object.
        temperature -- The temperature in C.
        dew_point -- The dew point in C.
        pressure -- the pressure in hPa
        """
        site = site.lower()
        assert isinstance(valid_time, datetime)

        if temperature is None or dew_point is None:
            return
        
        if dew_point > temperature:
            print(f"Humidity > 100: {site} {valid_time} {temperature}C / {dew_point}C")

        if valid_time.minute <= 30:
            vtime = valid_time - timedelta(minutes=valid_time.minute)
        else:
            vtime = valid_time + timedelta(minutes=60 - valid_time.minute)

        self.db_conn.execute(
            "INSERT OR REPLACE INTO obs (site_id, valid_time, temperature, dew, pres) VALUES (?, ?, ?, ?, ?)",
                (site, vtime, temperature, dew_point, pressure))
        return

    def get_hourly_data_for(self, site, starting=None, ending=None):
        """Get hourly temperature, dew point, and pressure data in a time range.

        # Arguments
        site - the site id, e.g. 'kmso', 'kgeg'
        starting - datetime.datetime for the start. If this is None it will grab the earliest 
            available data. If it has a value it will only retrieve data after this date (if it 
            exists in the database).
        ending - datetime.datetime for the end. If this is None it will grab up to the latest 
            available data. If it has a value it will only retrieve data before this date (if it 
            exists in the database).


        # Returns a generator object that yields tuples of 
        (site, valid time, temperature in C, dew point in C, pressure in hPa)
        """

        if starting is not None:
            assert isinstance(starting, datetime)
            starting = starting.strftime("'%Y-%m-%d %H:%M:%S'")

        if ending is not None:
            assert isinstance(ending, datetime)
            ending = ending.strftime("'%Y-%m-%d %H:%M:%S'")

        assert site is not None
        site = site.lower()

        query = """
            SELECT
                site_id,
                datetime(valid_time) as vt,
                temperature,
                dew,
                pres
            FROM obs
            """

        if site is not None:
            query += "           WHERE site_id = '%s' \n" % site

        if starting is not None:
            query += "                AND valid_time >= %s\n" % starting

        if ending is not None:
            query += "                AND valid_time <= %s\n" % ending

        
        query += "    ORDER BY vt ASC"

        c = self.db_conn.cursor()
        c.execute(query)

        vals = (list(r) for r in c.fetchall())

        def string_to_datetime(x_):
            long_date = datetime.strptime(x_, "%Y-%m-%d %H:%M:%S")
            return long_date

        def map_to_data(row):
            site = row[0]
            vt = string_to_datetime(row[1])
            temperature = row[2]
            dew = row[3]
            pres = row[4]

            if site is not None and vt is not None and temperature is not None and dew is not None:
                # Note that pres may be None at this point.
                return (site, vt, temperature, dew, pres)
            else:
                return None

        vals = (map_to_data(x) for x in vals)
        vals = (x for x in vals if x is not None)
        return vals


if __name__ == "__main__":

    import os
    import shutil

    # Test
    test_root = os.path.dirname(os.path.realpath(__file__))
    test_root = os.path.join(test_root, "test_workspace")

    # Remove the path from previous tests.
    if os.path.exists(test_root):
        print("Removing old test data...")
        shutil.rmtree(test_root)

    with Archive(root=test_root) as arch:
        print("First Instance\n")

    with Archive(root=test_root) as arch:
        print("Second Instance\n")
        arch.add_update_observation("kmso", datetime(2020, 7, 15, 0, 25), 70, 8, 800)
        arch.add_update_observation("kmso", datetime(2020, 7, 15, 1, 25), 70, 15, 800)
        arch.add_update_observation("kmso", datetime(2020, 7, 15, 2, 25), 70, 20)
        arch.add_update_observation("kmso", datetime(2020, 7, 15, 3, 25), 70, 25)
        arch.add_update_observation("kmso", datetime(2020, 7, 15, 4, 25), 70, 30)
        arch.add_update_observation("kmso", datetime(2020, 7, 15, 5, 25), 70, 34)
        arch.add_update_observation("kmso", datetime(2020, 7, 15, 6, 25), 70, 38)
        arch.add_update_observation("kmso", datetime(2020, 7, 15, 7, 25), 70, 40)
        arch.add_update_observation("kmso", datetime(2020, 7, 15, 8, 25), 70, 42)
        arch.add_update_observation("kmso", datetime(2020, 7, 15, 9, 25), 70, 44)
        arch.add_update_observation("kmso", datetime(2020, 7, 15, 10, 25), 70, 48)
        arch.add_update_observation("kmso", datetime(2020, 7, 15, 11, 25), 70, 50)
        arch.add_update_observation("kmso", datetime(2020, 7, 15, 12, 25), 70, 50)
        arch.add_update_observation("kmso", datetime(2020, 7, 15, 13, 25), 70, 54)
        arch.add_update_observation("kmso", datetime(2020, 7, 15, 14, 25), 70, 48)
        arch.add_update_observation("kmso", datetime(2020, 7, 15, 15, 25), 70, 40)
        arch.add_update_observation("kmso", datetime(2020, 7, 15, 16, 25), 70, 35)
        arch.add_update_observation("kmso", datetime(2020, 7, 15, 17, 25), 70, 30)
        arch.add_update_observation("kmso", datetime(2020, 7, 15, 18, 25), 70, 25)
        arch.add_update_observation("kmso", datetime(2020, 7, 15, 19, 25), 70, 20)
        arch.add_update_observation("kmso", datetime(2020, 7, 15, 20, 25), 70, 18)
        arch.add_update_observation("kmso", datetime(2020, 7, 15, 21, 25), 70, 14)
        arch.add_update_observation("kmso", datetime(2020, 7, 15, 22, 25), 70, 11)
        arch.add_update_observation("kmso", datetime(2020, 7, 15, 23, 25), 70, 9)
        arch.add_update_observation("kmso", datetime(2020, 7, 16, 0, 25), 70, 8)
        arch.add_update_observation("kmso", datetime(2020, 7, 16, 1, 25), 70, 15)
        arch.add_update_observation("kmso", datetime(2020, 7, 16, 2, 25), 70, 20)
        arch.add_update_observation("kmso", datetime(2020, 7, 16, 3, 25), 70, 25)
        arch.add_update_observation("kmso", datetime(2020, 7, 16, 4, 25), 70, 30)
        arch.add_update_observation("kmso", datetime(2020, 7, 16, 5, 25), 70, 34)
        arch.add_update_observation("kmso", datetime(2020, 7, 16, 6, 25), 70, 38)
        arch.add_update_observation("kmso", datetime(2020, 7, 16, 7, 25), 70, 40)
        arch.add_update_observation("kmso", datetime(2020, 7, 16, 8, 25), 70, 42)
        arch.add_update_observation("kmso", datetime(2020, 7, 16, 9, 25), 70, 44)
        arch.add_update_observation("kmso", datetime(2020, 7, 16, 10, 25), 70, 48)
        arch.add_update_observation("kmso", datetime(2020, 7, 16, 11, 25), 70, 50)
        arch.add_update_observation("kmso", datetime(2020, 7, 16, 12, 25), 70, 50)
        arch.add_update_observation("kmso", datetime(2020, 7, 16, 13, 25), 70, 54)
        arch.add_update_observation("kmso", datetime(2020, 7, 16, 14, 25), 70, 48)
        arch.add_update_observation("kmso", datetime(2020, 7, 16, 15, 25), 70, 40)
        arch.add_update_observation("kmso", datetime(2020, 7, 16, 16, 25), 70, 35)
        arch.add_update_observation("kmso", datetime(2020, 7, 16, 17, 25), 70, 30)
        arch.add_update_observation("kmso", datetime(2020, 7, 16, 18, 25), 70, 25)
        arch.add_update_observation("kmso", datetime(2020, 7, 16, 19, 25), 70, 20)
        arch.add_update_observation("kmso", datetime(2020, 7, 16, 20, 25), 70, 18)
        arch.add_update_observation("kmso", datetime(2020, 7, 16, 21, 25), 70, 14)
        arch.add_update_observation("kmso", datetime(2020, 7, 16, 22, 25), 70, 11)
        arch.add_update_observation("kmso", datetime(2020, 7, 16, 23, 25), 70, 9)

        vals = arch.get_hourly_data_for('kmso')
        for r in vals:
            print(r)

