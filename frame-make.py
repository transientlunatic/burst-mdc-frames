import sys,os,glob
import os.path
import otter

import pandas

#--------------------------------------#
# Startup information.
# ------
# This should be moved to a
# configuration file post haste
#--------------------------------------#

ifos = ['H1', 'L1']
runs = ['O1']
calibration = 'C01'

inj_families = ['ga'] # ['sg','ga'] # The types of injections for which frames must be produced.

xml_folder = '/home/daniel.williams/data/mdc/O1/xml'    # The location for the XML files
mdc_folder = '/home/daniel.williams/data/mdc/O1/frames' # The place where MDC frames will be stored once they've been generated.

report_location = '/home/daniel.williams/public_html/reports/mdc/o1'


# Look-up table to convert from the short to the long names for each injection family.
inj_families_names = {'ga' : 'Gaussian',
                      'sg' : 'SineGaussian',
                      'wnb': 'BTLWNB'
                      }





import numpy

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt

plt.style.use('ggplot')

from glue.ligolw import ligolw, utils, lsctables
lsctables.use_in(ligolw.LIGOLWContentHandler)

import lalburst, lalsimulation, lalmetaio

from pycbc.types import TimeSeries


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


# First we need to form a list of the data frames for which we need the MDC frames.


frames = pandas.DataFrame(columns=('run', 'ifo', 'start time', 'duration'))

i = 0
for run in runs:
    for ifo in ifos:
        print run, ifo
        root = path = "/archive/frames/{}/hoft_{}/{}".format(run, calibration, ifo)
        for r,d,f in os.walk(path):
            for file in f:
                file = file[:-4]   # strip the file extension out
                sttime = file[14:24]
                dur = file[25:]
                i+= 1
                frames.loc[i] = [run, ifo, sttime, dur]
                #report.write_row('Frame: st. {} dur. {}'.format(sttime, dur))

#report.write_row(df.to_html(classes="table"))



for family in inj_families:
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

        # Write a link into the main report for this injection family.
        report.write_row("<a href={}>{}</a>".format(subreport_loc, injection))
        
        # Load each injection xml file. Each starts with the family short name.
        # Get table
        xmldoc = utils.load_filename(injection, contenthandler=ligolw.LIGOLWContentHandler)
        sim_burst_tbl = lsctables.SimBurstTable.get_table(xmldoc)
        
        for sim in sim_burst_tbl:
            sim.time_geocent_gps += 1e-9*sim.time_geocent_gps_ns
        sims = filter(lambda s: s.waveform == inj_families_names[family], sim_burst_tbl)

        # Start the histograms of each of the file's parameters.
        report_inj.write_header(2,"Parameter histograms")
        fig, ax = plt.subplots(1,len(prm_hist_type[inj_families_names[family]]), sharey=True)   
        
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
            ax[i].semilogx()
            #ax[i].set_ylim([5e2, 1e5]);
            i+=1
        # Add the file to the report.
        report_inj.write_plot(figure=fig)
        plt.close()


        # Attempt to write out the frame files.
        # This requires access to lalsuite and to the 2014 Fall review branch
        # to run e.g.
        # lalapps_simburst_to_frame --simburst-file sg_f100_q9_linear_rescaled.xml --ifos [L1,H1] --gps-start 1125878400 --duration 4095

        for frame_info in frames:
            frame_loc = mdc_folder+"/"+family+"/"+injection[:-7]+"/"+frame_info['start time'][:5]+"/"+injection[:-7]+"-"+frame_info['start time']+"-"+frame_info['duration']
            # First check if the frame has already been made.
            if not os.path.isfile(frameloc):
                report_inj.write_row(frameloc)
            
        
        report_inj.write_footer()


report.write_footer()
