import sys,os,glob
import otter

import pandas

ifos = ['H1', 'L1']
runs = ['O1']
calibration = 'C01'

inj_families = ['ga'] # ['sg','ga'] # The types of injections for which frames must be produced.
inj_families_names = ['Gaussian'] # ['SineGaussian', 'Gaussian']

xml_folder = '/home/daniel.williams/data/mdc/O1/xml'    # The location for the XML files
mdc_folder = '/home/daniel.williams/data/mdc/O1/frames' # The place where MDC frames will be stored once they've been generated.

import numpy

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt

plt.style.use('ggplot')

from glue.ligolw import ligolw, utils, lsctables
lsctables.use_in(ligolw.LIGOLWContentHandler)

import lalburst, lalsimulation, lalmetaio

from pycbc.types import TimeSeries


report = otter.Otter('/home/daniel.williams/public_html/frametest.html', 
	{'title':"Test Report",
	'author':'Daniel Williams',
	'subtitle':''})


prm_hist_type = {
    "SineGaussian": { "hrss": "log",
            "time_geocent_gps": None,
            #"pol_ellipse_e": None,
            #"pol_ellipse_angle": None,
            "psi": None,
            "ra": None,
            "dec": None,
            #"frequency": "log",                                                                                                                                                                               
            #"q": "log" 
                      },                                                                                                                                                                                                            "Gaussian": { "hrss": "log",                                                                                                                                                                                                                  "time_geocent_gps": None,                                                                                                                                                                                                       #"pol_ellipse_e": None,                                                                                                                                                                                                         #"pol_ellipse_angle": None,    
                  "psi": None,
                  "ra": None,
                                                                                                                                                                                                                                                  "dec": None,
                                                                                                                                                                                                                                                  #"frequency": "log",                                                                                                                                                                                                 
                                                                                                                                                                                                                                                  #"q": "log"                                                                                                                                                                                                          
                      },  
    "BTLWNB": { "hrss": "log",
             "time_geocent_gps": None,
             "ra": None,
             "dec": None,             #"bandwidth": None,                                                                                                                                                                                         
             #"frequency": None,                                                                                                                                                                                      
             #"duration": None
                }             
    }



report.write_page_header()
#report.write_row("Hello world.")

for family in zip(inj_families, inj_families_names):
    os.chdir(xml_folder)
    for injection in glob.glob("{}*".format(family[0])):
        # Get table                                                                                                                                                                                                             
        xmldoc = utils.load_filename(injection, contenthandler=ligolw.LIGOLWContentHandler)
        sim_burst_tbl = lsctables.SimBurstTable.get_table(xmldoc)
        
        for sim in sim_burst_tbl:
            sim.time_geocent_gps += 1e-9*sim.time_geocent_gps_ns
        sims = filter(lambda s: s.waveform == family[1], sim_burst_tbl)

        fig, ax = plt.subplots(1,len(prm_hist_type[family[1]]), sharey=True)   
        
        i = 0
        for prm, htype in prm_hist_type[family[1]].iteritems():
            plt.setp(ax[i].get_xticklabels(), visible=True, rotation='vertical');

            fig.suptitle(injection)

            # Load the numbers in
            prms = [getattr(s, prm) for s in sims]

            ax[i].set_title(prm)
            if htype == "log":
                logbins = (min(prms), max(prms))
                logbins = numpy.logspace(numpy.log10(logbins[0]), numpy.log10(logbins[1]), 100, base=10)
                
                ax[i].hist(prms, bins=logbins, log=True, histtype='stepfilled')
                ax[i].semilogx()
                ax[i].set_ylim([5e2, 1e5]);
            else:
                ax[i].hist(prms, bins=100, histtype='stepfilled', log=True)
                ax[i].set_ylim([5e2, 1e5]);

            i+=1
        #fig.tight_layout()
        report.write_plot(figure=fig)
        plt.close()


df = pandas.DataFrame(columns=('run', 'ifo', 'start time', 'duration'))

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
                df.loc[i] = [run, ifo, sttime, dur]
                #report.write_row('Frame: st. {} dur. {}'.format(sttime, dur))

#report.write_row(df.to_html(classes="table"))

report.write_footer()
