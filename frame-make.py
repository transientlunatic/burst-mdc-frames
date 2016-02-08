import sys,os,glob
import os.path
import otter
import re
import pandas


import lalburst, lalsimulation, lalmetaio
from pylal.antenna import response
from logmake import *

#--------------------------------------#
# Startup information.
# ------
# This should be moved to a
# configuration file post haste
#--------------------------------------#

ifos = ['H1', 'L1']
runs = ['O1']
calibration = 'C01'

lasttime = 1136762880
#lasttime = 1137213440

inj_families = ['ga']# 'sg'] # ['sg','ga'] # The types of injections for which frames must be produced.

xml_folder = '/home/daniel.williams/data/mdc/O1/xml'    # The location for the XML files
mdc_folder = '/home/daniel.williams/data/mdc/O1/frames' # The place where MDC frames will be stored once they've been generated.

report_location = '/home/daniel.williams/public_html/reports/mdc/o1'


# Look-up table to convert from the short to the long names for each injection family.
inj_families_names = {'ga' : 'Gaussian',
                      'sg' : 'SineGaussian',
                      'wnb': 'BTLWNB'
                      }


# GravEN Header
header = "#\tGravEn_SimID\tSimHrss\tSimEgwR2\tGravEn_Ampl\tInternal_x\tInternal_phi\tExternal_x\tExternal_phi\tExternal_psi\tFrameGPS\tEarthCtrGPS\tSimName\tSimHpHp\tSimHcHc\tSimHpHc\tGEO\tGEOctrGPS\tGEOfPlus\tGEOfCross\tGEOtErr\tGEOcErr\tGEOsnr\tH1\tH1ctrGPS\tH1fPlus\tH1fCross\tH1tErr\tH1cErr\tH1snr\tH2\tH2ctrGPS\tH2fPlus\tH2fCross\tH2tErr\tH2cErr\tH2snr\tL1\tL1ctrGPS\tL1fPlus\tL1fCross\tL1tErr\tL1cErr\tL1snr\tV1\tV1ctrGPS\tV1fPlus\tV1fCross\tV1tErr\tV1cErr\tV1snr"


import numpy

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
font = {'size'   : 1}
matplotlib.rc('font', **font)

plt.style.use('ggplot')

from glue.ligolw import ligolw, utils, lsctables
lsctables.use_in(ligolw.LIGOLWContentHandler)

import lalburst, lalsimulation, lalmetaio

from pycbc.types import TimeSeries

def mkdir(path):
    sub_path = os.path.dirname(path)
    if not os.path.exists(sub_path):
        mkdir(sub_path)
    if not os.path.exists(path):
        os.mkdir(path)


report = otter.Otter(report_location+"/index.html", 
	{'title':"O1 Burst MDC",
	'author':'Daniel Williams',
	'subtitle':''})

# Lists of the quantities which will be plotted in the histograms for each injection family.
prm_hist_type = {
    "SineGaussian": { "hrss": "log",
            "time_geocent_gps": None,
            "psi": None,
            "ra": None,
            "dec": None,
    },
    "Gaussian": { "hrss": "log",
                  "time_geocent_gps": None,
                  "psi": None,
                  "ra": None,
                  "dec": None,
    },  
    "BTLWNB": { "hrss": "log",
                "time_geocent_gps": None,
                "ra": None,
                "dec": None,
    }             
}


report.write_page_header()
report.write_breadcrumb([('Burst MDCs', '../'),('O1', '#')])


def write_burst_mdc_row( row, start=0 ):
    """
    Fill in a template row of a BurstMDC style (GravEn) log.

    Template:
    template = simulation_id, hrss, egw, graven_amp, internal_x, internal_phi, external_x, external_phi, external_psi, frame_gps, earth_center_gps, waveform, simhphp, simhxhx, simhphx, g1fp, g1fx, 0, 0, 0, h1fp, h1fx, 0, 0, 0, h2fp, h2fx, 0, 0, 0, l1fp, l1fx, 0, 0, 0, v1fp, v1fx, 0, 0, 0
    """
    
    sim_id = row.simulation_id

    print row.simulation_id
    
    sim_hrss = row
    # FIXME: What are these?
    graven_amp = 0
    internal_x = 0
    internal_phi = 0
    ###########
    external_x = row.dec
    external_phi = row.ra
    external_psi = row.psi
    frame_gps = start
    earth_center_gps = float(row.time_geocent_gps)
    sim_name = row.waveform
    # The factor of 1.8597e-21 is to convert from units of M_s/pc^2 to SI units
    egw = (row.egw_over_rsquared or 0)*1.8597e-21
    # FIXME: This is wrong by the light travel time from earth center to det
    g1fp, g1fx, _, _ = response(earth_center_gps, row.ra, row.dec, 0, row.psi, 'radians', "G1")
    h1fp, h1fx, _, _ = response(earth_center_gps, row.ra, row.dec, 0, row.psi, 'radians', "H1")
    h2fp, h2fx, _, _ = response(earth_center_gps, row.ra, row.dec, 0, row.psi, 'radians', "H2")
    l1fp, l1fx, _, _ = response(earth_center_gps, row.ra, row.dec, 0, row.psi, 'radians', "L1")
    v1fp, v1fx, _, _ = response(earth_center_gps, row.ra, row.dec, 0, row.psi, 'radians', "V1")
    _, simhphp, simhxhx, simhphx = measure_hrss(row)
    hrss = row.hrss or 0
    template = "%d\t%g\t%g\t%f\t%f\t%f\t%f\t%f\t%f\t%10.10f\t%10.10f\t%s\t%g\t%g\t%g\tGEO\t%f\t%f\t%f\t%f\t%f\tH1\t%f\t%f\t%f\t%f\t%f\tH2\t%f\t%f\t%f\t%f\t%f\tL1\t%f\t%f\t%f\t%f\t%f\tV1\t%f\t%f\t%f\t%f\t%f" % (int(row.simulation_id), hrss, egw, graven_amp, internal_x, internal_phi, external_x, external_phi, external_psi, frame_gps, earth_center_gps, row.waveform, simhphp, simhxhx, simhphx, g1fp, g1fx, 0, 0, 0, h1fp, h1fx, 0, 0, 0, h2fp, h2fx, 0, 0, 0, l1fp, l1fx, 0, 0, 0, v1fp, v1fx, 0, 0, 0)
    return template

def write_burst_mdc_log(fname, rows):
    """
    Write out a set of BurstMDC (GravEn) rows into fname.
    """
    f = open(fname, "w")
    print >>f, header
    for row in rows:
        print >>f, row
    f.close()

def make_logfile(frameloc, frame_info, sim_burst_tbl):
    if not os.path.isfile(frameloc+'.gwf.log'):
        frame_info = frames.loc[j]
        start, end = float(frame_info['start time']), float(frame_info['start time'])+float(frame_info['duration'])
        mdc_log = []
        for row in sim_burst_tbl:
            #print row.time_geocent_gps, start, end
            if row.time_geocent_gps > end: break
            elif row.time_geocent_gps < start: continue
            mdc_log.append(write_burst_mdc_row(row, start))
        write_burst_mdc_log(frameloc+'.gwf.log', mdc_log)


    
# First we need to form a list of the data frames for which we need the MDC frames.

frames = pandas.DataFrame(columns=('run', 'ifo', 'start time', 'duration'))

i = 0
for run in runs:
    for ifo in ifos:
        print run, ifo
        root = path = "/archive/frames/{}/hoft/{}".format(run, ifo)
        for r,d,f in os.walk(path):
            for file in f:
                file = file[:-4]   # strip the file extension out
                sttime = file[14:24]
                dur = file[25:]
                # Check if the frame already exists for another ifo, and add this to the list of ifos instead
                coinc = frames[(frames['start time']==sttime) &  (frames['duration']==dur)]
                if len(coinc)>0:
                    ind=coinc.index.values[0]
                    frames.ix[ind]['ifo'] = numpy.append(numpy.array(coinc['ifo'])[0], numpy.array(ifo))
                else:
                    frames.loc[i] = [run, numpy.array(ifo), sttime, dur]
                    i+= 1
#report.write_row('Frame: st. {} dur. {}'.format(sttime, dur))

#report.write_row(df.to_html(classes="table"))



for family in inj_families:
    report.write_header(2,inj_families_names[family])
    report._write("<a href='{}'></a>".format(family))
    # Work on each family of injections sequentially.
    os.chdir(xml_folder)
    
    for injection in glob.glob("{}*".format(family)):
        # Each injection gets its own report
        subreport_loc="{}/{}.html".format(family,injection)
        report_inj = otter.Otter(report_location+"/"+subreport_loc, 
	                         {'title':"O1 Burst MDC",
	                          'author':'Daniel Williams',
	                          'subtitle':'{}'.format(injection)}
        )
        report_inj.write_page_header()
        report_inj.write_breadcrumb([('Burst MDCs', '../../'),('O1','../'),(inj_families_names[family], '../index.html#{}'.format(family)),(injection,'#')])

        
        # Write a link into the main report for this injection family.
        report.write_row("<a href={}>{}</a>".format(subreport_loc, injection))

        
        # Load each injection xml file. Each starts with the family short name.
        # Get table
        #xmldoc = utils.load_filename(injection, contenthandler=ligolw.LIGOLWContentHandler)
        #sim_burst_tbl = lsctables.SimBurstTable.get_table(xmldoc)

        sim_burst_table = lalburst.SimBurstTableFromLIGOLw(injection, None, None)
        waveforms = []
        while True:
            waveforms.append(sim_burst_table)
            if sim_burst_table.next is None: break
            sim_burst_table = sim_burst_table.next
        
        #for sim in waveforms:
        #    sim.time_geocent_gps += 1e-9*sim.time_geocent_gps_ns
        sims = filter(lambda s: s.waveform == inj_families_names[family], waveforms)

        # Start the histograms of each of the file's parameters.
        report_inj.write_header(2,"Parameter histograms")
        fig, ax = plt.subplots(1,len(prm_hist_type[inj_families_names[family]]), sharey=True, figsize=(12,3))   
        
        i = 0
        for prm, htype in prm_hist_type[inj_families_names[family]].iteritems():
            plt.setp(ax[i].get_xticklabels(), visible=True, rotation='vertical');
            fig.suptitle(injection)
            # Load the numbers in
            prms = [getattr(s, prm) for s in sims]
            ax[i].set_title(prm)
            logbins = (min(prms), max(prms))
            #logbins = numpy.logspace(numpy.log10(logbins[0]), numpy.log10(logbins[1]), 100, base=10)
            ax[i].hist(prms, bins=100, log=True, histtype='stepfilled')
            #ax[i].semilogx()
            ax[i].set_ylim([1e2, 5e5]);
            i+=1
        plt.tight_layout()
        # Add the file to the report.
        report_inj.write_plot(figure=fig)
        plt.close()

        report_inj.write_warning('info', 'These injections can be found in {}'.format(mdc_folder+"/"+family+"/"+injection[:-16]+"/"))



        # Attempt to write out the frame files.
        # This requires access to lalsuite and to the 2014 Fall review branch
        # to run e.g.
        # lalapps_simburst_to_frame --simburst-file sg_f100_q9_linear_rescaled.xml --ifos [L1,H1] --gps-start 1125878400 --duration 4095

        for j in xrange(len(frames)):
            frame_info = frames.loc[j]
            if int(frame_info['duration']) < 20: continue
            if int(frame_info['start time']) > lasttime: continue
            
            ifosstr =  "".join(re.findall("[a-zA-Z]+", str(frame_info['ifo'])))
            
            frameloc = mdc_folder+"/"+family+"/"+injection[:-16]+"/"+frame_info['start time'][:5]+"/"+ifosstr+"-"+inj_families_names[family]+"-"+frame_info['start time']+"-"+frame_info['duration']

            
            # First check if the frame has already been made.
            if not os.path.isfile(frameloc+".gwf"):
                
                path = mdc_folder+"/"+family+"/"+injection[:-16]+"/"+frame_info['start time'][:5]
                mkdir(path)
                
                # Move into this folder so that we can produce the frame here
                os.chdir(path)

                try: 
                    os.system("/home/daniel.williams/repositories/lalsuite/lalapps/src/power/lalapps_simburst_to_frame --simburst-file {} --ifos [{}] --gps-start {} --duration {}".format(xml_folder+"/"+injection, 
                                                                                                                             str(frame_info['ifo']).replace("'","").replace(" ",",").replace("[","").replace("]",""),
                                                                                                                             frame_info['start time'], 
                                                                                                                             frame_info['duration']))
                except:
                    if not os.path.isfile(frameloc):
                        report_inj.write_warning('danger', "Failed frame creation: {}".format(frame_info['start time']))

            # Make the GravEn log file for the frame
            make_logfile(frameloc, frame_info, sim_burst_tbl)

                        
        os.chdir(xml_folder)
        report_inj.write_footer()


report.write_footer()
