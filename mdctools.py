from glue.ligolw import ligolw, utils, lsctables
lsctables.use_in(ligolw.LIGOLWContentHandler);
import numpy
import lalburst, lalsimulation, lalmetaio
from pylal.antenna import response

from pylal.date import XLALTimeDelayFromEarthCenter
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
from pylal import inject 
import numpy as np
import pandas as pd


class MDCSet():
    def __init__(self, detectors, simtable, name=None):
        """
        Represents an MDC set, stored in an XML SimBurstTable file.
        
        Parameters
        ----------
        detectors : list 
            A list of detector names where the injections should be made
        simtable : str
            The filepath to a simbursttable xml file.
        """
        self.detectors = detectors
        sim_burst_table = lalburst.SimBurstTableFromLIGOLw(simtable, None, None)
        self.waveforms = []
        self.strains = []
        self.times = []
        while True:
            self.waveforms.append(sim_burst_table)
            self.strains.append(self._generate_burst(sim_burst_table))
            self.times.append(sim_burst_table.time_geocent_gps)
            if sim_burst_table.next is None: break
            sim_burst_table = sim_burst_table.next
        if name: 
            self.name = name
        else:
            self.name = self._simID(0)
            
        self.times = np.array(self.times)
        
    def _generate_burst(self, row,rate=16384.0):
        """
        Generate the burst described in a given row, so that it can be 
        measured.
        
        Parameters
        ----------
        row : int
            The row number of the waveforms to be measured
        rate : float
            The sampling rate of the signal, in Hz. Defaults to 16384.0Hz
            
        Returns
        -------
        hp : 
            The strain in the + polarisation
        hx : 
            The strain in the x polarisation
        hp0 : 
            A copy of the strain in the + polarisation
        hx0 : 
            A copy of the strain in the x polarisation
        """
        #row = self.waveforms[row]
        self.swig_row = lalburst.CreateSimBurst()
        for a in lsctables.SimBurstTable.validcolumns.keys():
            try:
                setattr(self.swig_row, a, getattr( row, a ))
            except AttributeError: continue # we didn't define it
            except TypeError: 
                print a, getattr(row,a)
                continue # the structure is different than the TableRow
        hp, hx = lalburst.GenerateSimBurst(self.swig_row, 1.0/rate)
        # FIXME: Totally inefficent --- but can we deep copy a SWIG SimBurst?
        # DW: I tried that, and it doesn't seem to work :/
        hp0, hx0 = lalburst.GenerateSimBurst(self.swig_row, 1.0/rate)
        return hp, hx, hp0, hx0 
    
    def _getDetector(self, det):
        """
        A method to return a LALDetector object corresponding to a detector's
        X#-style name, e.g. 'H1' as the Hanford 4km detector.
        
        Parameters
        ----------
        det : str
            A string describing the detector in the format letter-number, e.g
            "H1" would be the Hanford 4km detector, "L1" would be the 
            Livingston 4km, and so-forth.
            
        Returns
        -------
        detector : LALDetector
            The LAL object describing the detector
        """
        # create detector-name map 
        detMap = {'H1': 'LHO_4k', 'H2': 'LHO_2k', 'L1': 'LLO_4k', 
                  'G1': 'GEO_600', 'V1': 'VIRGO', 'T1': 'TAMA_300'} 
        try: 
            detector=detMap[det] 
        except KeyError: 
             raise ValueError, "ERROR. Key %s is not a valid detector name." % (det) 

        # get detector 
        if detector not in inject.cached_detector.keys(): 
              raise ValueError, "%s is not a cached detector.  "\
                    "Cached detectors are: %s" % (det, inject.cached_detector.keys()) 
        return inject.cached_detector[detector]

    def _timeDelayFromGeocenter(self, detector, ra, dec, gpstime):
        """
        Calculate the time delay between the geocentre and a given detector
        for a signal from some sky location.
        
        Parameters
        ----------
        detector : str
            A string describing the detector, e.g. H1 is the Hanford 4km 
            detector.
        ra : float
            The right-ascension of the observation in radians
        dec : float
            The declination of the obser
        """
        if isinstance(detector, str): detector = self._getDetector(detector)
        gpstime = LIGOTimeGPS(float(gpstime))
        return XLALTimeDelayFromEarthCenter(detector.location, ra, dec, gpstime)
    
    def _simID(self, row):
        """
        Generate a name for an injection set in the format expected by cWB
        
        Parameters
        ----------
        row : SimBurst
            The simburst table row describing the injection
        """
        row = self.waveforms[row]
        name = ''

        inj_families_names = {'ga' : 'Gaussian',
                          'sg' : 'SineGaussian',
                          'wnb': 'BTLWNB'
                          }
        inj_families_abb = dict((v,k) for k,v in inj_families_names.iteritems())
        name += '{}{:.3f}'.format(inj_families_abb[row.waveform].upper(), row.duration*1e3).replace('.','d')

        return name
    
    def _measure_hrss(self, row, rate=16384.0):
        """
        Measure the various components of hrss (h+^2, hx^2, hphx) for a given 
        input row. This is accomplished by generating the burst and calling 
        the SWIG wrapped  XLALMeasureHrss in lalsimulation. 
        
        Parameters
        ----------
        row : int
            The row number of the waveforms to be measured
        rate : float
            The sampling rate of the signal, in Hz. Defaults to 16384.0Hz
            
        Returns
        -------
        hrss : float
            The measured hrss of the waveform amplitude: sqrt(|Hp|^2 + |Hx|^2)
        hphp : float
            The hrss of the + polarisation only.
        hxhx : float
            The hrss of the x polarisation only.
        hphx : float
            The hrss of |HpHx| 
        """
        hp, hx, hp0, hx0 = self.strains[row]
        hp0.data.data *= 0
        hx0.data.data *= 0

        # H+ hrss only
        hphp = lalsimulation.MeasureHrss(hp, hx0)**2
        # Hx hrss only
        hxhx = lalsimulation.MeasureHrss(hp0, hx)**2
        # sqrt(|Hp|^2 + |Hx|^2)
        hrss = lalsimulation.MeasureHrss(hp, hx)

        hp.data.data = numpy.abs(hx.data.data) + numpy.abs(hp.data.data)
        # |H+Hx|
        hphx = (lalsimulation.MeasureHrss(hp, hx0)**2 - hrss**2)/2
        return hrss, hphp, hxhx, hphx
    
    def _measure_egw_rsq(self, row, rate=16384.0):
        """
        Measure the energy emitted in gravitational waves divided 
        by the distance squared in M_solar / pc^2. This is accomplished 
        by generating the burst and calling the SWIG wrapped  
        XLALMeasureHrss in lalsimulation. 
        
        Parameters
        ----------
        row : int
            The row number of the waveforms to be measured
        rate : float
            The sampling rate of the signal, in Hz. Defaults to 16384.0Hz
            
        Returns 
        -------
        egw : float
            The energy emitted in gravitational waves divided 
            by the distance squared in M_solar / pc^2.
        """
        hp, hx, _, _ = self.strains[row]
        return lalsimulation.MeasureEoverRsquared(hp, hx)
    
    def _responses(self, row):
        """
        Calculate the antenna repsonses for each detector to the waveform.
        
        Parameters
        ----------
        row : int
            The row number of the waveforms to be measured
            
        Returns
        -------
        responses : list of lists
            A list containing the lists of antenna responses, with the first 
            element of each list containing the detector acronym.
        """
        output = []
        row = self.waveforms[row]
        for detector in self.detectors:
            time = row.time_geocent_gps + self._timeDelayFromGeocenter(detector, row.ra, row.dec, row.time_geocent_gps)
            time = np.float64(time)
            rs = response(time, row.ra, row.dec, 0, row.psi, 'radians', detector)
            output.append([detector, time, rs[0], rs[1]]   )
        return output
    
    def gravEn_row(self, row, frame):
        """
        Produces a gravEn-style log row for a row of the simBurstTable.
        
        Parameters
        ----------
        row : int
            The row number of the waveforms to be measured
            
        Returns
        -------
        str
            A string in the gravEn format which describes the injection.
        """
        strains = self._measure_hrss(row)
        rowname = self._simID(row)
        responses = self._responses(row)
        energy = self._measure_egw_rsq(row)
        row = self.waveforms[row]
        output = [] 
        output.append(self.name)                  # GravEn_SimID
        output.append(strains[0])                 # SimHrss
        output.append(energy)                     # SimEgwR2
        output.append(strains[0])                 # GravEn_Ampl
        output.append(0)                          # Internal_x (currently not implemented)
        output.append(0)                          # Intenal_phi ('')
        output.append(row.ra)                     # External_x
        output.append(row.dec)                    # External_phi
        output.append(row.psi)                    # External_psi
        output.append(frame.start)                # FrameGPS
        output.append(row.time_geocent_gps)       # EarthCtrGPS
        output.append(rowname)                    # SimName
        output.append(strains[1])                 # SimHpHp
        output.append(strains[2])                 # SimHcHc
        output.append(strains[3])                 # SimHpHp
        output.append(" ".join(" ".join(map(str,l)) for l in responses))
        return ' '.join(str(e) for e in output)

class Frame():
    """
    Represents a frame, in order to prepare the injection frames
    """
    def __init__(self, start, duration, ifo):
        self.start = start
        self.duration = duration
        self.end = self.start + duration
        self.ifos = ifo
        
    def __repr__(self):
        out = ''
        out += "MDC Frame \n"
        for ifo in self.ifos:
            out += "{} {} {} \n".format(ifo, self.start, self.duration)
        return out
    
    def get_rowlist(self,mdcs):
        """
        Return the rows from an MDCs which correspond to this frame.
        
        Parameters
        ----------
        mdcs : MDCSet object
            The set of MDCs from which the rows are to be found.
        """
        return np.where((mdcs.times<self.end)&(mdcs.times>self.start))[0]
    
    def calculate_n_injections(self, mdcs):
        return len(mdcs.times[(mdcs.times<self.end)&(mdcs.times>self.start)])
    
    def generate_log(self,mdc):
        log = '#  GravEn_SimID  SimHrss  SimEgwR2  GravEn_Ampl  Internal_x  Internal_phi  External_x  External_phi External_psi  FrameGPS  EarthCtrGPS  SimName  SimHpHp  SimHcHc  SimHpHc  H1       H1ctrGPS        H1fPlus        H1fCross    L1       L1ctrGPS        L1fPlus        L1fCross\n'
        rowlist = self.get_rowlist(mdc)
        for row in rowlist:
            log += mdc.gravEn_row(row, self)
            log += "\n"
        return log


class FrameSet():

    def __init__(self, frame_list):
        """
        A collection of frames.

        Parameters
        ----------
        frame_list : str
            The filespath of a CSV file containing the list of frames, 
            and the parameters required to produce them: the start and 
            duration times, and the interferometers they describe.
        """

        self.frames = []
        frame_list = pd.read_csv(frame_list)
        for frame in frame_list.iterrows():
            frame = frame[1]
            ifos = frame['ifo'].replace("['",'').replace("']",'').replace("'",'').split(' ')
            frame = Frame(frame['start time'],frame['duration'],ifos)
            self.frames.append(frame)
        
    def full_logfile(self, mdc, location):
        """
        Produce a log file for the entire frame set 
        """
        full_log = ''
        for frame in self.frames:
            full_log += frame.generate_log(mdc)
            
        text_file = open(location, "w")
        text_file.write(full_log)
        text_file.close()

