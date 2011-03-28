''' Create a database connection and retrieve the time sampling information 
    The constructor connects to the default database and returns a cursor to the database.
    The getTimeMagSQL takes a filter, ra/dec pair, and a boolean to which determines whether
    to query the original cronos.92 pointings or the dithered ones.
    Created September 24 2007 by K. Simon Krughoff University of Washington.
    Modified:
    March 2009 by K. Simon Krughoff krughoff@astro.washington.edu
'''
import MySQLdb
import numpy as num
import math
class DB:
    ''' Class for handling database connections '''
    host = "lsst-db.astro.washington.edu"
    user = "lsst"
    passwd = "lsst"
    dbase = "lsst_pointings"
    def __init__(self):
        self.__enter__()

    def __enter__(self):
        self.db = MySQLdb.connect(host=self.host, user=self.user, passwd=self.passwd, db=self.dbase)
        self.cursor = self.db.cursor()
        return self

    def __exit__(self, *dumArgs):
        if self.db:
            self.cursor.close()
            self.db.close()
            self.db = None
            self.cursor = None

    def close(self):
        '''Close the db connection; may be safely called even if already closed'''
        self.__exit__()

    def isOpen(self):
        return self.db != None
   
    def getTimeMagSQL(self, ra, dec, filt, doDith=False, version="opsim3_61"):
        ''' Retrieve timing from default database for a single ra/dec pair in a particular filter from
            any of the OpSim runs.
        Inputs:
        - ra: a sequence of RA positions in degrees
        - dec: a sequence of declination positions in degrees
        - filt: a string describing the filter for the sampling.  Currently the sloan u, g, r, i, z, and y filters are acceptable values
        - doDith: boolean to state which database to use, the un-dithered (False) or dithered pointings (True)
        - version: Version of the operation simulator to use.  Currently, there are three: Cronos92 (old), OpSim5_72 and OpSim1_29 (new).
        Return:
        - time: a sequence of simulated observation in MJD
        - m5: a sequence of the 5 sigma limiting magnitude associated with each observation
        '''
        lversion = version.lower()
        m5col = ""
        simtab = ""
        if(lversion == "cronos92"):
            m5col = "m5"
            simtab = "pointing_meta"
        elif(lversion == "opsim5_72"):
            m5col = "5sigma_ps"
            simtab = "output_opsim5_72"
        elif(lversion == "opsim1_29"):
            m5col = "5sigma_ps"
            simtab = "output_opsim1_29"
        elif(lversion == "opsim3_61"):
            m5col = "5sigma_ps"
            simtab = "output_opsim3_61"
        else:
            m5col = "5sigma_ps"
            simtab = "output_opsim3_61"
            print "Didn't recognize the run name (%s).  Setting to latest: OpSim3_61"%version
           
        if ra < 0 or ra > 360.:
            ra = ra%360.
        if dec > 90:
           dec = 90.
        if dec < -90:
           dec = -90.
        time = []
        #List to contain the time sampling
        m5 = []
        #List to contain the 5 sigma limiting magnitudes for the pointings
        queryparts = []
        pointing_radius_deg = 1.75
        deg2rad = math.pi/180.0
        rad2deg = 180.0/math.pi
        #Constant to convert degrees to radians
        pointing_radius_rad = pointing_radius_deg*deg2rad

        piover2 = math.pi/2.0
        #draw box around the pointing
        decmin = dec - pointing_radius_deg;
        decmax = dec + pointing_radius_deg;
        if decmin < -90:
            decmin = -90.
        if decmax > 90:
            decmax = 90.
        declimstr = "between %f and %f"%(decmin*deg2rad,decmax*deg2rad)

        decorr = decmin
        if math.fabs(decmin) > math.fabs(decmax):
            decorr = decmin
        else:
            decorr = decmax
        ralimstr = []
        if math.fabs(decorr) >= 90.:
            ralimstr = ["between %f and %f"%(0.*deg2rad, 360.*deg2rad)]
        else:
            ramin = ra - pointing_radius_deg/math.cos(decorr*deg2rad)
            ramax = ra + pointing_radius_deg/math.cos(decorr*deg2rad) 
            if ramin < 0 and ramax > 360:
                ralimstr = [ "between %f and %f"%(0.*deg2rad, 360.*deg2rad) ]
            elif ramin < 0:
                ralimstr.append("between %f and %f"%((360.+ ramin)*deg2rad, 360.*deg2rad))
                ralimstr.append("between %f and %f"%(0.*deg2rad, ramax*deg2rad))
            elif ramax > 360:
                ralimstr.append("between %f and %f"%(0.*deg2rad, (ramax - 360)*deg2rad))
                ralimstr.append("between %f and %f"%(ramin*deg2rad, 360.*deg2rad))
            else:
                ralimstr = [ "between %f and %f"%(ramin*deg2rad, ramax*deg2rad) ]

        if doDith:
            ralimstr = " or hexdithra ".join(ralimstr)
            queryparts.append("select distinct(b.expMJD), b.%s from (select a.* from "%(m5col))
            queryparts.append("(select hexdithra, hexdithdec from %s where hexdithdec %s and (hexdithra %s) group by hexdithra, hexdithdec) a where "%(simtab, declimstr, ralimstr))
            queryparts.append("acos(")
            queryparts.append("sin(%f - hexdithdec)*cos(hexdithra)*sin(%f - %f*(%f))*cos(%f*%f) + "%(piover2, piover2, deg2rad, dec, deg2rad, ra))
            queryparts.append("sin(%f - hexdithdec)*sin(hexdithra)*sin(%f - %f*(%f))*sin(%f*%f) + "%(piover2, piover2, deg2rad, dec, deg2rad, ra))
            queryparts.append("cos(%f - hexdithdec)*cos(%f - %f*%f) "%(piover2, piover2, deg2rad, dec))
            queryparts.append(")< %f) a, %s b where "%(pointing_radius_rad, simtab))
            queryparts.append("a.hexdithra = b.hexdithra and a.hexdithdec = b.hexdithdec and filter = \'%s\' order by b.expMJD"%(filt))
            #Query string to do the 3 space dot product of the pointing specified by ra/dec and all nondithered field centers.
            #Selects all field centers within 0.03054/deg2rad = 1.74981 degrees of the ra/dec pair
        else:
            ralimstr = " or fieldra ".join(ralimstr)
            queryparts.append("select distinct(b.expMJD), b.%s from (select a.* from "%(m5col))
            queryparts.append("(select fieldra, fielddec from %s where fielddec %s and (fieldra %s) group by fieldra, fielddec) a where "%(simtab, declimstr, ralimstr))
            queryparts.append("acos(")
            queryparts.append("sin(%f - fielddec)*cos(fieldra)*sin(%f - %f*(%f))*cos(%f*%f) + "%(piover2, piover2, deg2rad, dec, deg2rad, ra))
            queryparts.append("sin(%f - fielddec)*sin(fieldra)*sin(%f - %f*(%f))*sin(%f*%f) + "%(piover2, piover2, deg2rad, dec, deg2rad, ra))
            queryparts.append("cos(%f - fielddec)*cos(%f - %f*%f) "%(piover2, piover2, deg2rad, dec))
            queryparts.append(")< %f) a, %s b where "%(pointing_radius_rad, simtab))
            queryparts.append("a.fieldra = b.fieldra and a.fielddec = b.fielddec and filter = \'%s\' order by b.expMJD"%(filt))
            #Query string to do the 3 space dot product of the pointing specified by ra/dec and all dithered field centers.
            #Selects all field centers within 0.03054/deg2rad = 1.74981 degrees of the ra/dec pair
        query = "".join(queryparts)
        #Join pieces of the query string
        self.cursor.execute(query)
        result = self.cursor.fetchall()
        #Execute query and return all results
        #Check if query returns.  If not, return empty arrays.
        if len(result) == 0:
            time = num.asarray(time)
            m5 = num.asarray(m5)
        else:
            #Read results into arrays.  The star is very important.
            time, m5 = zip(*result)
            time = num.asarray(time)
            m5 = num.asarray(m5)
        return time, m5
