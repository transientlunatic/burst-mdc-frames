from glue.ligolw import ligolw, utils, lsctables
lsctables.use_in(ligolw.LIGOLWContentHandler);

import lalburst, lalsimulation, lalmetaio
from pylal.antenna import response

injection = 'ga_d0100_rescaled.xml.gz'
start = 1126621184

xmldoc = utils.load_filename(injection, contenthandler=ligolw.LIGOLWContentHandler)
sim_burst_tbl = lsctables.SimBurstTable.get_table(xmldoc)


def write_burst_mdc_log(fname, rows):
    """
    Write out a set of BurstMDC (GravEn) rows into fname.
    """
    f = open(fname, "w")
    print >>f, header
    for row in rows:
        print >>f, row
    f.close()

def measure_egw_rsq(row, rate=16384.0):
    """
    Measure the energy emitted in gravitational waves divided by the distance squared in M_solar / pc^2. This is accomplished by generating the burst and calling the SWIG wrapped  XLALMeasureHrss in lalsimulation. Thus, the row object should be a SWIG wrapped SimBurst object. Rate is the sampling rate in Hz (default is 16kHz).
    """
    swig_row = lalburst.CreateSimBurst()
    for a in lsctables.SimBurstTable.validcolumns.keys():
        try:
            setattr(swig_row, a, getattr( row, a ))
        except AttributeError: continue # we didn't define it
        except TypeError: continue # the structure is different than the TableRow
    hp, hx = lalburst.GenerateSimBurst(swig_row, 1.0/rate)
    return lalsimulation.MeasureEoverRsquared(hp, hx)

def measure_hrss(row, rate=16384.0):
    """
    Measure the various components of hrss (h+^2, hx^2, hphx) for a given input row. This is accomplished by generating the burst and calling the SWIG wrapped  XLALMeasureHrss in lalsimulation. Thus, the row object should be a SWIG wrapped SimBurst object. Rate is the sampling rate in Hz (default is 16kHz).
    """
    swig_row = lalburst.CreateSimBurst()
    for a in lsctables.SimBurstTable.validcolumns.keys():
        try:
            setattr(swig_row, a, getattr( row, a ))
        except AttributeError: continue # we didn't define it
        except TypeError: continue # the structure is different than the TableRow
    hp, hx = lalburst.GenerateSimBurst(swig_row, 1.0/rate)
    # FIXME: Totally inefficent --- but can we deep copy a SWIG SimBurst?
    hp0, hx0 = lalburst.GenerateSimBurst(swig_row, 1.0/rate)
    hp0.data.data *= 0
    hx0.data.data *= 0

    # H+ hrss only
    hphp = lalsimulation.MeasureHrss(hp, hx0)**2
    # Hx hrss only
    hxhx = lalsimulation.MeasureHrss(hp0, hx)**2
    # sqrt(|Hx|^2 + |Hx|^2)
    hrss = lalsimulation.MeasureHrss(hp, hx)

    hp.data.data = numpy.abs(hx.data.data) + numpy.abs(hp.data.data)
    # |H+Hx|
    hphx = (lalsimulation.MeasureHrss(hp, hx0)**2 - hrss**2)/2
    return hrss, hphp, hxhx, hphx


def write_burst_mdc_row( row, start=0, end=0):
    """
    Fill in a template row of a BurstMDC style (GravEn) log.

    Template:
    template = simulation_id, hrss, egw, graven_amp, internal_x, internal_phi, external_x, external_phi, external_psi, frame_gps, earth_center_gps, waveform, simhphp, simhxhx, simhphx, g1fp, g1fx, 0, 0, 0, h1fp, h1fx, 0, 0, 0, h2fp, h2fx, 0, 0, 0, l1fp, l1fx, 0, 0, 0, v1fp, v1fx, 0, 0, 0
    """
    
    sim_id = row.simulation_id
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
    earth_center_gps = float(row.get_time_geocent())
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

mdc_log = []
for row in sim_burst_tbl:
    if row.time_geocent_gps > end: break
    elif row.time_geocent_gps < start: continue
    mdc_log.append(write_burst_mdc_row(row, start))
write_burst_mdc_log('injections.log', mdc_log)
