import mdctools

injection="ga_d0100_rescaled.xml.gz"

mdcs = mdctools.MDCSet(['H1', 'L1'], injection)

o1 = mdctools.FrameSet('frame_list.dat')
o1.full_logfile(mdcs, 'logfile.txt')
