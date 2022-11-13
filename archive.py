"""An archive of point observations.

 This file is a library to serve as an interface to an sqlite3 file storing observation data.
"""
from collections import deque, namedtuple
from collections.abc import Iterable
from datetime import datetime, date, timedelta
import math as m
from pathlib import Path
import sqlite3
import tqdm


from watermark import watermark
print(f"Import versions for libraries imported by {__name__}")
print(watermark(iversions=True, globals_=globals()))

Location = namedtuple('Location', ['loc_id', 'stn_num', 'site_id', 'skip_import', 'latitude',
    'longitude', 'elevation', 'count', 'description'])

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
        #print("Archive root is", root)

        # Check if the databases exists.
        db_file = root / DB_FILE
        #if not db_file.exists():
        #    print("No database found, one will be created.")

        self.db_conn = sqlite3.connect(db_file,
                                       detect_types=sqlite3.PARSE_DECLTYPES
                                       | sqlite3.PARSE_COLNAMES)
        #print("Connected to", db_file)

        # Set up the tables.
        c = self.db_conn.cursor()
        c.executescript("""
                -- Hourly Data, no averaging
                CREATE TABLE IF NOT EXISTS obs (
                    loc_id      INTEGER   NOT NULL,
                    valid_time  TIMESTAMP NOT NULL,
                    temp        REAL      NOT NULL, -- Temperature in C
                    dew         REAL      NOT NULL, -- Dew point in C
                    vpd         REAL      NOT NULL, -- Vapor Pressure Deficit in hPa
                    qc_level    INTEGER   NOT NULL, -- Lower is better.
                    UNIQUE(loc_id, valid_time)
                );

                CREATE INDEX IF NOT EXISTS faster_queries ON obs (loc_id, valid_time);

                -- Locations and sites.
                CREATE TABLE IF NOT EXISTS locations (
                    loc_id      INTEGER PRIMARY KEY,
                    stn_num     INTEGER NOT NULL,
                    site_id     TEXT    NOT NULL,
                    description TEXT,
                    skip_import INTEGER DEFAULT 0, -- Default to importing data.
                    latitude    REAL    NOT NULL,  -- Latitude in degrees.
                    longitude   REAL    NOT NULL,  -- Longitude in degrees.
                    elevation   REAL    NOT NULL,  -- Elevation in meters.
                    count       INTEGER NOT NULL   -- Number of observations in the database.
                    );

                CREATE INDEX IF NOT EXISTS faster_joins ON locations (site_id, loc_id);
            """)

        # Track if there is modifications.
        self._modified = False

        # Load all the locations.
        self.__locations = self.__load_locations()
        self.__next_loc_id = None

        # New stations to skip
        self.__new_skips = []

        return

    def __load_locations(self):
        '''Get a list of all the locations in the database.

        # Returns
        A dictionary with stn_num as keys and a list of namedtuple objects of type Location for 
        values.
        '''
        c = self.db_conn.cursor()
        query = """
            SELECT 
                loc_id,
                stn_num, 
                site_id, 
                skip_import, 
                latitude, 
                longitude, 
                elevation, 
                count,
                description
            FROM locations;"""

        def map_row_to_object(row):
            loc_id = int(row[0])
            stn_num = int(row[1])
            site_id = row[2]
            skip_import = row[3] > 0
            latitude = float(row[4])
            longitude = float(row[5])
            elevation = float(row[6])
            count = int(row[7])
            description = row[8]

            return Location(
                    loc_id=loc_id, 
                    stn_num=stn_num, 
                    site_id=site_id, 
                    skip_import=skip_import, 
                    latitude=latitude, 
                    longitude=longitude, 
                    elevation=elevation, 
                    count=count,
                    description=description
                )

        d = {}
        for row in c.execute(query):
            loc = map_row_to_object(row)
            if loc.stn_num not in d:
                d[loc.stn_num] = []
            d[loc.stn_num].append(loc)

        c.close()
        del c
        return d
            

    def __query_count(self, loc_id):
        c = self.db_conn.cursor()
        c.execute("SELECT COUNT(loc_id) FROM obs WHERE loc_id=?", (loc_id,))
        count = c.fetchone()
        if count is None:
            return 0
        else:
            return count[0]


    def __store_locations(self):
        '''Save self.__locations in the database.'''
        locations = [loc for stn_num, loc_list in self.__locations.items() for loc in loc_list]

        for i in range(len(locations)):
            count = self.__query_count(locations[i].loc_id)
            locations[i] = locations[i]._replace(count=count)

        query = """
            INSERT INTO locations (
                loc_id,
                stn_num, 
                site_id, 
                skip_import, 
                latitude, 
                longitude, 
                elevation, 
                count,
                description
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT DO UPDATE SET count = excluded.count
        """

        with self.db_conn:
            c = self.db_conn.cursor()
            c.executemany(query, locations)
            c.close()

        return


    def __get_next_loc_id(self):
        '''Get the next value for loc_id.'''
        if self.__next_loc_id is None:
            if len(self.__locations) == 0:
                return 0
            max_id = max(loc.loc_id for key, val in self.__locations.items() for loc in val)
            self.__next_loc_id = max_id + 1

        next_id = self.__next_loc_id
        self.__next_loc_id += 1
        return next_id

    def __add_location(self, loc):
        self._modified = True

        if loc.stn_num in self.__locations:
            for old_loc in self.__locations[loc.stn_num]:
                assert loc.loc_id != old_loc
            self.__locations[loc.stn_num].append(loc)
        else:
            self.__locations[loc.stn_num] = [loc]

        return


    def __get_location_for(self, stn_num, latitude, longitude, elevation, site_id, description):
        '''Retrieve a Location object matching the above criteria, or create a new one.'''
        skip = stn_num in self.__new_skips

        if stn_num in self.__locations:
            for i in range(len(self.__locations[stn_num])):
                c= self.__locations[stn_num][i]
                if c.latitude == latitude and c.longitude == longitude and c.elevation == elevation:
                    self._modified = True

                    self.__locations[stn_num][i] = self.__locations[stn_num][i]._replace(skip_import=skip)
                    if site_id is not None:
                        self.__locations[stn_num][i] = self.__locations[stn_num][i]._replace(site_id=site_id)

                    if description is not None:
                        self.__locations[stn_num][i] = self.__locations[stn_num][i]._replace(description=description)

                    return self.__locations[stn_num][i]

        if site_id is None:
            site_id = "missing"

        new_location = Location(self.__get_next_loc_id(), stn_num, site_id, skip, latitude,
                longitude, elevation, 0, description)

        self.__add_location(new_location)

        return new_location


    def __enter__(self):
        return self

    def __exit__(self, *args):
        if self._modified:
            self.__store_locations()

        self.db_conn.close()

        return

    def __get_xhour_table_names(self, hours):
        '''Get the table names and index names for a moving average that is hours long.'''
        hours = int(hours)

        assert (hours % 10) == 0 and hours > 0

        data_table_name = f"obs{hours}"
        index_name = f"faster_queries_{hours}"

        return (data_table_name, index_name)

    def add_skip_stations(self, stn_nums):
        '''Make sure these stations are in the skip list so no data for them is ever added.'''
        if isinstance(stn_nums, Iterable):
            stn_nums = list(stn_num for stn_num in stn_nums)
        else:
            stn_nums = [stn_nums]

        self.__new_skips.extend(stn_nums)

        return

    def add_all(self, data_gen):
        """Add all the observations in data_gen to the archive.

        Arguments:
        data_gen -- an iterable that returns tuples of (station number, station id, coordinates,
            valid time, temperature, dew point, vapor pressure deficit, QC code).
        """
        self._modified = True

        # Get the location info
        matched = (
                (
                    self.__get_location_for(stn_num, latitude, longitude, elevation, site_id, name), 
                    vt,
                    temp,
                    dew,
                    vpd,
                    qc
                ) for 
                    stn_num, 
                    site_id,
                    (latitude, longitude, elevation),
                    vt, 
                    temp, 
                    dew, 
                    vpd, 
                    qc,
                    name
                in data_gen)

        matched = ((loc.loc_id, vt, temp, dp, vpd, qc) for 
                loc, vt, temp, dp, vpd, qc in matched if not loc.skip_import)

        with self.db_conn:

            self.db_conn.executemany(
                """
                INSERT OR IGNORE INTO obs (
                    loc_id,
                    valid_time, 
                    temp,
                    dew,
                    vpd,
                    qc_level
                ) VALUES (?, ?, ?, ?, ?, ?)
                """, 
                matched
            )

        self.__store_locations()

        return


    def get_all_sites(self):
        """Get a list of all the sites in the archive."""
        sites = set()
        for stn_num, locations in self.__locations.items():
            for location in locations:
                if not location.skip_import and location.site_id != 'missing':
                    sites.add(location.site_id)

        return sorted(sites)


    def remove_stations(self, stn_nums):
        '''Remove all entries for a station number (NOT site id).'''
        self.__store_locations()
        self._modified = True

        if isinstance(stn_nums, Iterable):
            stn_nums = tuple(stn_num for stn_num in stn_nums)
        else:
            stn_nums = (stn_nums,)

        stn_nums = tuple(stn_num for stn_num in stn_nums if stn_num in self.__locations)

        loc_ids = tuple((loc.loc_id,)
                for stn_num in stn_nums 
                for loc in self.__locations[stn_num]
                )

        num_rows = 0
        with self.db_conn:
            num_rows += self.db_conn.executemany("DELETE FROM obs WHERE loc_id=?", loc_ids).rowcount
            num_rows += self.db_conn.executemany(
                    "UPDATE locations SET skip_import = 1 WHERE loc_id=?", loc_ids).rowcount

        self.__locations = self.__load_locations()

        return num_rows


    def __get_location_info_for_site(self, site):
        '''Get the site_id, average latitude, longitude, and elevation for a site.'''

        count = 0
        sum_lat = 0.0
        sum_lon = 0.0
        sum_elev = 0.0

        for _, locations in self.__locations.items():
            for location in locations:
                if location.site_id == site:
                    count += 1
                    sum_lat += location.latitude
                    sum_lon += location .longitude
                    sum_elev += location.elevation

        if count == 0:
            return None
        return (sum_lat / count, sum_lon / count, sum_elev / count)


    def get_hourly_vpd_data_for(self, site, starting=None, ending=None):
        """Get hourly vapor pressure deficit data in a time range.

        # Arguments
        site - the site id, e.g. 'kmso', 'kgeg'
        starting - datetime.datetime for the start. If this is None it will grab the earliest 
            available data. If it has a value it will only retrieve data after this date (if it 
            exists in the database).
        ending - datetime.datetime for the end. If this is None it will grab up to the latest 
            available data. If it has a value it will only retrieve data before this date (if it 
            exists in the database).


        # Returns a generator object that yields tuples of 
        (site, valid time, vapor pressure deficit in hPa, latitude, longitude, elevation, 
            number of matching obs)
        """

        if starting is not None:
            assert isinstance(starting, datetime)
            starting = starting.strftime("'%Y-%m-%d %H:%M:%S'")
            start_str = f"                AND valid_time >= {starting}"
        else:
            start_str = ""

        if ending is not None:
            assert isinstance(ending, datetime)
            ending = ending.strftime("'%Y-%m-%d %H:%M:%S'")
            end_str = f"                AND valid_time <= {ending}"
        else:
            end_str = ""

        assert site is not None
        site = site.lower()

        avg_coords = self.__get_location_info_for_site(site)
        assert avg_coords is not None
        latitude, longitude, elevation = avg_coords

        c = self.db_conn.cursor()
        c.executescript(
                f"""
                    DROP TABLE IF EXISTS hourly;

                    CREATE TEMP TABLE hourly as 
                    SELECT
                        valid_time, 
                        vpd, 
                        qc_level 
                    FROM obs JOIN locations on obs.loc_id = locations.loc_id 
                    WHERE site_id = '{site}'
                        {start_str}
                        {end_str}
                    ORDER BY valid_time ASC;
                """)


        query = f"""
            SELECT
                ob1.valid_time, 
                AVG(vpd), 
                COUNT(*)
            FROM hourly as ob1
            JOIN (
                SELECT valid_time, MIN(qc_level) as qc
                FROM hourly
                GROUP BY valid_time) as ob2
            ON ob1.valid_time = ob2.valid_time 
                AND ob2.qc = ob1.qc_level
            GROUP BY ob1.valid_time
            ORDER BY ob1.valid_time ASC
            """

        c.execute(query)

        def map_to_data(row):
            vt = datetime.strptime(row[0],"%Y-%m-%d %H:%M:%S")
            vpd = row[1]
            count = row[2]

            if vt is not None and vpd is not None:
                return (site, vt, vpd, latitude, longitude, elevation, count)
            else:
                return None

        vals = (map_to_data(x) for x in c)
        vals = (x for x in vals if x is not None)
        return vals

    def build_xhour_vpd_averages(self):
        '''Build tables with running averages of VPD.'''

        periods = (10, 100, 1000, 10000)

        for hours in periods:
            table_name, index_name = self.__get_xhour_table_names(hours)

            self.db_conn.executescript(
                f"""
                    DROP TABLE IF EXISTS {table_name};

                    CREATE TABLE IF NOT EXISTS {table_name} (
                        site_id     TEXT      NOT NULL, -- e.g. kmso, kgpi
                        valid_time  TIMESTAMP NOT NULL,
                        vpd         REAL      NOT NULL  -- Vapor Pressure Deficit in hPa
                    );

                    CREATE INDEX IF NOT EXISTS {index_name} ON {table_name} (site_id, valid_time);
                """)

        min_counts = tuple(int(0.9 * hours) for hours in periods)

        def avg_deque(dq_vpd, min_count):
            count = len(dq_vpd)
            if count < min_count:
                return None

            avg_vpd = sum(dq_vpd) / count

            return avg_vpd

        def build_data_gen(site):

            buffers_vt = tuple(deque([]) for _ in periods)
            buffers_vpd = tuple(deque([]) for _ in periods)

            data_gen = self.get_hourly_vpd_data_for(site)

            for _, vt, vpd, _, _, _, _ in data_gen:

                # Add new data to the buffer
                for buffer_vt, buffer_vpd in zip(buffers_vt, buffers_vpd):
                    buffer_vt.append(vt)
                    buffer_vpd.append(vpd)

                # Remove old data
                earliests = tuple(vt - timedelta(hours=hours) for hours in periods)
                for earliest, buffer_vt, buffer_vpd in zip(earliests, buffers_vt, buffers_vpd):

                    while len(buffer_vt) > 0 and buffer_vt[0] < earliest:
                        buffer_vt.popleft()
                        buffer_vpd.popleft()

                avg_vpds = tuple(avg_deque(buffer_vpd, min_count)
                        for buffer_vpd, min_count in zip(buffers_vpd, min_counts))

                yield (site, vt, avg_vpds)

            return
                    
        print("Building moving average tables.")

        table_names = (self.__get_xhour_table_names(hours) for hours in periods)
        table_names = tuple(x[0] for x in table_names)

        queries = tuple(
                f""" 
                    INSERT INTO {table_name} (site_id, valid_time, vpd) VALUES (?, ?, ?)
                """ for table_name in table_names)

        sites = self.get_all_sites()
        pbar = tqdm.tqdm(sites)
        for site in pbar:

            pbar.set_description(site)
            
            for site, vt, avg_vpds in build_data_gen(site):
                for query, avg_vpd in zip(queries, avg_vpds):
                    if avg_vpd is not None:
                        self.db_conn.execute(query, (site, vt, avg_vpd))
            self.db_conn.commit()

        return

    def get_xhourly_average_vpd_data_for(self, hours, site, starting=None, ending=None):
        """Get the averaged hourly vapor pressure deficit data in a time range.

        # Arguments
        hours - how long of an average. 
        site - the site id, e.g. 'kmso', 'kgeg'
        starting - datetime.datetime for the start. If this is None it will grab the earliest 
            available data. If it has a value it will only retrieve data after this date (if it 
            exists in the database).
        ending - datetime.datetime for the end. If this is None it will grab up to the latest 
            available data. If it has a value it will only retrieve data before this date (if it 
            exists in the database).


        # Returns a generator object that yields tuples of 
        (site, valid time, average vapor pressure deficit in hPa, average latitude, 
        average longitude, average elevation)
        """

        if starting is not None:
            assert isinstance(starting, datetime)
            starting = starting.strftime("'%Y-%m-%d %H:%M:%S'")
            start_str = f"                AND valid_time >= {starting}"
        else:
            start_str = ""

        if ending is not None:
            assert isinstance(ending, datetime)
            ending = ending.strftime("'%Y-%m-%d %H:%M:%S'")
            end_str = f"                AND valid_time <= {ending}"
        else:
            end_str = ""

        assert site is not None
        site = site.lower()

        avg_coords = self.__get_location_info_for_site(site)
        assert avg_coords is not None
        latitude, longitude, elevation = avg_coords

        table_name, _ = self.__get_xhour_table_names(hours)

        query = f"""
            SELECT
                valid_time, 
                vpd 
            FROM {table_name}
            WHERE site_id = '{site}'
                {start_str}
                {end_str}
            ORDER BY valid_time ASC
            """

        c = self.db_conn.cursor()
        c.execute(query)
        vals = (list(r) for r in c)

        def map_to_data(row):
            vt = row[0]
            vpd = row[1]

            if vt is not None and vpd is not None:
                return (site, vt, vpd, latitude, longitude, elevation)
            else:
                return None

        vals = (map_to_data(x) for x in vals)
        vals = (x for x in vals if x is not None)
        return vals

    def output_sites_kml(self, kml_file, log_file):
        '''Create a kml file with info about the sites.'''

        kml_file = open(kml_file, "wt")

        #
        # Initialize the KML output file.
        #
        kml_file.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
        kml_file.write("<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n")
        kml_file.write("<Document>\n")

        kml_file.write("<Style id=\"circle\">\n")
        kml_file.write("<PolyStyle>\n<fill>0</fill>\n<outline>1</outline>\n</PolyStyle>\n")
        kml_file.write("<LineStyle>\n<color>ff0000cc</color>\n")
        kml_file.write("<colorMode>normal</colorMode>\n<width>2</width>\n</LineStyle>\n")
        kml_file.write("<IconStyle>\n<scale>0</scale>\n</IconStyle>\n")
        kml_file.write("</Style>\n")
    
        kml_file.write("<Style id=\"circleballoon\">\n")
        kml_file.write("<IconStyle>\n<Icon>\n")
        kml_file.write("<href>http://maps.google.com/mapfiles/kml/paddle/grn-blank.png</href>\n")
        kml_file.write("</Icon>\n</IconStyle>\n")
        kml_file.write("<LabelStyle><scale>0</scale></LabelStyle>\n")
        kml_file.write("</Style>\n")
    
        kml_file.write("<Style id=\"missing\">\n")
        kml_file.write("<IconStyle>\n<Icon>\n")
        kml_file.write("<href>http://maps.google.com/mapfiles/kml/paddle/red-stars.png</href>\n")
        kml_file.write("</Icon>\n</IconStyle>\n")
        kml_file.write("<LabelStyle><scale>0</scale></LabelStyle>\n")
        kml_file.write("</Style>\n")

        kml_file.write("<Style id=\"skipped\">\n")
        kml_file.write("<IconStyle>\n<Icon>\n")
        kml_file.write("<href>http://maps.google.com/mapfiles/kml/shapes/forbidden.png</href>\n")
        kml_file.write("</Icon>\n</IconStyle>\n")
        kml_file.write("<LabelStyle><scale>0</scale></LabelStyle>\n")
        kml_file.write("</Style>\n")

        #
        # Organize the stations by their type
        #
        missing_stations = []
        skip_stations = []
        stations_by_id = {}
        for _, locations_list in self.__locations.items():
            for loc in locations_list:
                if loc.skip_import:
                    skip_stations.append(loc)
                elif loc.site_id == 'missing':
                    missing_stations.append(loc)
                else:
                    if loc.site_id not in stations_by_id:
                        stations_by_id[loc.site_id] = []
                    stations_by_id[loc.site_id].append(loc)

        #
        # Output the stations and their circles.
        #
        distance_warnings = []
        elevation_warnings = []
        kml_file.write("<Folder><open>1</open><visibility>1</visibility><name>Stations</name>\n");
        for site_id in sorted(stations_by_id.keys()):
            log_file.write(f"\n{site_id}\n")
            kml_file.write("<Folder><open>0</open><visibility>1</visibility>\n");
            kml_file.write(f"<name>{site_id}</name>\n")
    
            all_coords = []
            for loc in stations_by_id[site_id]:
                log_file.write(f"\t{loc.stn_num:12d} | {loc.count:8d}\n")
    
                all_coords.append((loc.latitude, loc.longitude, loc.elevation))

                log_file.write(f"\t\t\t{loc.latitude:8.5f} | {loc.longitude:10.5f} | {loc.elevation: 7.1f}\n")
                kml_file.write(f"<Placemark><name>{loc.stn_num}</name>\n")
                kml_file.write(f"<description>{loc.stn_num} - {loc.description}</description>\n")
                kml_file.write("<styleUrl>#circleballoon</styleUrl>\n<Point><coordinates>")
                kml_file.write(f"{loc.longitude},{loc.latitude},0</coordinates></Point></Placemark>\n")
    
            # The maximum distances between sets of coordinates.
            max_dist = 0.0
            max_elev_diff = 0.0
            farthest_pair = None
            for i, (lat, lon, elev) in enumerate(all_coords):
                if i < (len(all_coords) - 1):
                    for j in range(i + 1, len(all_coords), 1):
                        lat2, lon2, elev2 = all_coords[j]
                        dist = haversine_distance(lat, lon, lat2, lon2)

                        if max_dist < dist:
                            farthest_pair = ((lat, lon),(lat2, lon2))

                        elev_diff = abs(elev - elev2)
                        max_elev_diff = max(max_elev_diff, elev_diff)
                        max_dist = max(max_dist, dist)

            if farthest_pair is not None:
                (lat1, lon1), (lat2, lon2) = farthest_pair
                latc = (lat1 + lat2) / 2
                lonc = (lon1 + lon2) / 2
                circle_points = make_haversine_circle(latc, lonc, max_dist / 2)

                kml_file.write("<Placemark>\n")
                kml_file.write(f"<name>{site_id.upper()}</name>\n")
                kml_file.write("<styleUrl>#circle</styleUrl>\n")
                kml_file.write("<MultiGeometry>\n")
                kml_file.write(f"<Point><coordinates>{lonc},{latc},0</coordinates></Point>\n")
                kml_file.write("<Polygon>\n<outerBoundaryIs>\n")
                kml_file.write("<LinearRing>\n<coordinates>\n")
                for lat, lon in circle_points:
                    kml_file.write(f"{lon},{lat},0\n")
                kml_file.write("</coordinates>\n</LinearRing>\n")
                kml_file.write("</outerBoundaryIs>\n</Polygon>\n")
                kml_file.write("</MultiGeometry>\n</Placemark>\n")

            log_file.write(f"\t\t\tMax Distance is {max_dist:.3f}km\n")
            log_file.write(f"\t\t\tMax Elevation Difference is {max_elev_diff:.3f}m\n")
            log_file.write("\n\n")

            if max_dist > 1.0: # km
                distance_warnings.append((site_id, max_dist))

            if max_elev_diff > 3.0: # m
                elevation_warnings.append((site_id, max_elev_diff))

            kml_file.write("</Folder>\n")

        kml_file.write("</Folder>\n")

        distance_warnings.sort(key=lambda x: x[1], reverse=True)
        elevation_warnings.sort(key=lambda x: x[1], reverse=True)

        for site_id, max_dist in distance_warnings:
            log_file.write(f"WARNING: {site_id} has max distance of {max_dist:.3f}km between coordinates.\n")
        for site_id, max_elev_diff in elevation_warnings:
            log_file.write(
                f"WARNING: {site_id} has max elevation difference of {max_elev_diff:.3f}m between coordinates.\n")

        #
        # Output missing stations
        #
        kml_file.write("<Folder><open>1</open><visibility>1</visibility><name>Missing Stations</name>\n")
        for loc in missing_stations:
            log_file.write(f"Could not assign {loc.stn_num}: {loc.site_id} {loc.latitude} ")
            log_file.write(f"{loc.longitude} {loc.elevation} {loc.description}\n")

            kml_file.write(f"<Placemark><name>{loc.stn_num}</name>\n")
            kml_file.write(f"<description>{loc.site_id}-{loc.stn_num}-{loc.description}</description>\n")
            kml_file.write("<styleUrl>#missing</styleUrl>\n<Point><coordinates>")
            kml_file.write(f"{loc.longitude},{loc.latitude},0</coordinates></Point></Placemark>\n")
    
        kml_file.write("</Folder>\n")

        #
        # Output skip stations
        #
        kml_file.write("<Folder><open>0</open><visibility>1</visibility><name>Skipped Stations</name>\n")
        for loc in skip_stations:
            log_file.write(f"Skipped {loc.stn_num}: {loc.site_id} {loc.latitude} ")
            log_file.write(f"{loc.longitude} {loc.elevation} {loc.description}\n")

            kml_file.write(f"<Placemark><name>{loc.stn_num}</name>\n")
            kml_file.write(f"<description>{loc.site_id}-{loc.stn_num}-{loc.description}</description>\n")
            kml_file.write("<styleUrl>#skipped</styleUrl>\n<Point><coordinates>")
            kml_file.write(f"{loc.longitude},{loc.latitude},0</coordinates></Point></Placemark>\n")
    
        kml_file.write("</Folder>\n")

        #
        # Put closing tags and close kml file.
        #
        kml_file.write("    </Document>\n")
        kml_file.write("</kml>\n")
        kml_file.close()



###################################################################################################
# Calculate the approximate distance between two points on the Earth's surface.
###################################################################################################
#
# This is useful for a quick check that I haven't grouped the stations poorly. Stations do move, 
# but it shouldn't be too far.
def haversine_distance(lat1, lon1, lat2, lon2):
    R = 6371 # Radius of Earth in km
    dLat = m.radians(lat2 - lat1)
    dLon = m.radians(lon2 - lon1)
    a = m.sin(dLat / 2)**2 + \
            m.cos(m.radians(lat1)) * m.cos(m.radians(lat2)) * m.sin(dLon / 2)**2
    c = 2 * m.atan2(m.sqrt(a), m.sqrt(1 - a))
    d =  R * c; # Distance in km
    return d

###################################################################################################
# Calculate the points in a circle around a point with a given distance using Haversine formula.
###################################################################################################
def make_haversine_circle(origin_lat, origin_lon, distance_km):
    # Number of points in the circle
    NUM_POINTS = 16

    # Radius of Earth in km
    R = 6371 

    lat1 = m.radians(origin_lat)
    lon1 = m.radians(origin_lon)

    outputs = []
    for i in range(NUM_POINTS):
        bearing = m.radians(float(i) / NUM_POINTS * 360.0)

        lat2 = m.asin(m.sin(lat1) * m.cos(distance_km / R) + 
                m.cos(lat1) *m.sin(distance_km / R) * m.cos(bearing))
        lon2 = lon1 + m.atan2(m.sin(bearing) * m.sin(distance_km / R) * m.cos(lat1),
                m.cos(distance_km / R) - m.sin(lat1) * m.sin(lat2))

        lat2 = m.degrees(lat2)
        lon2 = m.degrees(lon2)

        outputs.append((lat2, lon2))
    outputs.append(outputs[0])

    return outputs


