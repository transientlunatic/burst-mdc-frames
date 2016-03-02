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

import mdctools

o1 = mdctools.FrameSet('frame_list.dat')



inj_families = ['sg'] # ['sg','ga'] # The types of injections for which frames must be produced.

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


    
# First we need to form a list of the data frames for which we need the MDC frames.


# i = 0
# for run in runs:
#     for ifo in ifos:
#         print run, ifo
#         root = path = "/archive/frames/{}/hoft/{}".format(run, ifo)
#         for r,d,f in os.walk(path):
#             for file in f:
#                 file = file[:-4]   # strip the file extension out
#                 sttime = file[14:24]
#                 dur = file[25:]
#                 # Check if the frame already exists for another ifo, and add this to the list of ifos instead
#                 coinc = frames[(frames['start time']==sttime) &  (frames['duration']==dur)]
#                 if len(coinc)>0:
#                     ind=coinc.index.values[0]
#                     frames.ix[ind]['ifo'] = numpy.append(numpy.array(coinc['ifo'])[0], numpy.array(ifo))
#                 else:
#                     frames.loc[i] = [run, numpy.array(ifo), sttime, dur]
#                     i+= 1

frames = o1.frame_list

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
        # print "Making figures..."

        # sim_burst_table = lalburst.SimBurstTableFromLIGOLw(injection, None, None)
        # waveforms = []
        # while True:
        #     waveforms.append(sim_burst_table)
        #     if sim_burst_table.next is None: break
        #     sim_burst_table = sim_burst_table.next
        
        # sims = filter(lambda s: s.waveform == inj_families_names[family], waveforms)

        # # Start the histograms of each of the file's parameters.
        # report_inj.write_header(2,"Parameter histograms")
        # fig, ax = plt.subplots(1,len(prm_hist_type[inj_families_names[family]]), sharey=True, figsize=(12,3))   
        
        # i = 0
        # for prm, htype in prm_hist_type[inj_families_names[family]].iteritems():
        #     plt.setp(ax[i].get_xticklabels(), visible=True, rotation='vertical');
        #     fig.suptitle(injection)
        #     # Load the numbers in
        #     prms = [getattr(s, prm) for s in sims]
        #     ax[i].set_title(prm)
        #     logbins = (min(prms), max(prms))
        #     #logbins = numpy.logspace(numpy.log10(logbins[0]), numpy.log10(logbins[1]), 100, base=10)
        #     ax[i].hist(prms, bins=100, log=True, histtype='stepfilled')
        #     #ax[i].semilogx()
        #     ax[i].set_ylim([1e2, 5e5]);
        #     i+=1
        # plt.tight_layout()
        # # Add the file to the report.
        # report_inj.write_plot(figure=fig)
        # plt.close()

        # report_inj.write_warning('info', 'These injections can be found in {}'.format(mdc_folder+"/"+family+"/"+injection[:-16]+"/"))



        # Attempt to write out the frame files.
        # This requires access to lalsuite and to the 2014 Fall review branch
        # to run e.g.
        # lalapps_simburst_to_frame --simburst-file sg_f100_q9_linear_rescaled.xml --ifos [L1,H1] --gps-start 1125878400 --duration 4095

        # for j in xrange(len(frames)):
        #     frame_info = frames.loc[j]
        #     if int(frame_info['duration']) < 20: continue
        #     if int(frame_info['start time']) > lasttime: continue
            
        #     ifosstr =  "".join(re.findall("[a-zA-Z]+", str(frame_info['ifo'])))
            
             
        #     head_date = str(frame_info['start time'])[:5]
        #     #print mdc_folder +"/"+ family+"/" + injection[:-16] + "/" + head_date + "/" + ifosstr + inj_families_names[family],
        #     #print str(frame_info['start time']) + frame_info['duration']
        #     frameloc = mdc_folder+"/"+family+"/"+injection[:-16]+"/"+head_date+"/"+ifosstr+"-"+inj_families_names[family]+"-"+str(frame_info['start time'])+"-"+str(frame_info['duration'])

            
        #     # First check if the frame has already been made.
        #     if not os.path.isfile(frameloc+".gwf"):
                
        #         path = mdc_folder+"/"+family+"/"+injection[:-16]+"/"+frame_info['start time'][:5]
        #         mkdir(path)
                
        #         # Move into this folder so that we can produce the frame here
        #         os.chdir(path)

        #         try: 
        #             os.system("/home/daniel.williams/repositories/lalsuite/lalapps/src/power/lalapps_simburst_to_frame --simburst-file {} --ifos [{}] --gps-start {} --duration {}".format(xml_folder+"/"+injection, 
        #                                                                                                                      str(frame_info['ifo']).replace("'","").replace(" ",",").replace("[","").replace("]",""),
        #                                                                                                                      frame_info['start time'], 
        #                                                                                                                      frame_info['duration']))
        #         except:
        #             if not os.path.isfile(frameloc):
        #                 report_inj.write_warning('danger', "Failed frame creation: {}".format(frame_info['start time']))


        # First check if thelog has been made
        log_filename = "LAL"+injection[:-16]+".log"
        log_filepath = mdc_folder+"/"+family+"/"+injection[:-16]+"/input/"

        if not os.path.isfile(log_filepath+log_filename):            
            # Make the GravEn log file for the frame
            mdcs = mdctools.MDCSet(['H1', 'L1'], injection)
            # Make a folder to put it in the place cWB like
            mkdir(mdc_folder+"/"+family+"/"+injection[:-16]+"/"+"input/")
            # Make the file name
            print "Making log file for {}".format(injection[:-16])
            print "\t in {}".format(log_filepath+log_filename)
            o1.full_logfile(mdcs, log_filepath+log_filename)
            os.system("sort -k 11 {} > {}".format(log_filepath+log_filename, log_filepath+log_filename[:-4]+"_gps_sorted.log"))
            del(mdcs)
                        
        os.chdir(xml_folder)
        report_inj.write_footer()
        

report.write_footer()
