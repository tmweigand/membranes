import mdsys as md
import mdsys_io as io
import timestep as st
import timeseries as se
import numpy as np
import argparse
import sys


parser = argparse.ArgumentParser(description='Pass desired folder to data analysis tools.')
parser.add_argument('-i', '--in_folder', help ='folder to perform data analysis on')
parser.add_argument('--old', action="store_true")
args = parser.parse_args()


print('Building time series analysis object.')
sys.stdout.flush()
series_analyzer = se.AnalyzeTimeSeries(args.in_folder)
print('Beginning time series analysis.')
sys.stdout.flush()
series_analyzer.write_gfgs(is_old_format=args.old)
