#	Copyright (C) 2013 Martin D. Hurst
#	This program is free software: you can redistribute it and/or modify it 
#	under the terms of the GNU General Public License as published by the 
#	Free Software Foundation, either version 3 of the License, or (at your 
#	option) any later version. This program is distributed in the hope that 
#	it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
#	warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License along 
#	with this program. If not, see http://www.gnu.org/licenses/.
#
#	Martin D. Hurst mhurst@bgs.ac.uk

def get_means(infile, meanfile):

    #import modules
    import sys, os, numpy as np

    #infile = sys.argv[1]
    #meanfile = sys.argv[2]


    if (os.path.exists(os.getcwd() + "/" + infile) == False):
        raise sys.exit("Cannot find file " + infile) 
    else:
        print "\tGetting mean data for hilltops in " + infile + "...",

    ## open file
    f = open(infile,'r')
    lines = f.readlines()
    no_lines = len(lines)

    # data variables
    X = np.zeros(no_lines)
    Y = np.zeros(no_lines)
    seg = np.zeros(no_lines)
    cht = np.zeros(no_lines)
    slope = np.zeros(no_lines)
    relief = np.zeros(no_lines)
    Lh = np.zeros(no_lines)

    # read in data
    for i in range (1,no_lines):
        line = lines[i].strip().split(" ")
        
        X[i-1] = (line[0])
        Y[i-1] = (line[1])
        seg[i-1] = int(line[2])
        cht[i-1] = float(line[3])
        slope[i-1] = float(line[4])
        relief[i-1] = float(line[5])
        Lh[i-1] = float(line[6])

    f.close()

    # sort by segments
    max_seg = int(np.max(seg))

    #mean variables
    cht_mean = np.zeros(max_seg)
    slope_mean = np.zeros(max_seg)
    relief_mean = np.zeros(max_seg)
    Lh_mean = np.zeros(max_seg)

    for j in range (0,max_seg+1):
        count = 0
        seg_temp = np.zeros(no_lines)
        cht_temp = np.zeros(no_lines)
        slope_temp = np.zeros(no_lines)
        relief_temp = np.zeros(no_lines)
        Lh_temp = np.zeros(no_lines)

        for i in range(0,no_lines):
            if seg[i] == j:
                cht_temp[count] = cht[i]
                slope_temp[count] = slope[i] 
                relief_temp[count] = relief[i]
                Lh_temp[count] = Lh[i]
                count = count+1

        cht_mean[j-1] = np.mean(cht_temp[0:count-1])
        slope_mean[j-1] = np.mean(slope_temp[0:count-1])
        relief_mean[j-1] = np.mean(relief_temp[0:count-1])
        Lh_mean[j-1] = np.mean(Lh_temp[0:count-1])

    seg = np.arange(0,max_seg)

    #clear up for nans and slope

##    count = 0
##    for i in range (0,max_seg):
##        if (np.isnan(slope_mean[i]) == 0 and np.isnan(cht_mean[i]) == 0 and slope_mean[i] > 0 and cht_mean[i] < 0):
##            seg_temp[count] = seg[i]
##            cht_temp[count] = cht_mean[i]
##            slope_temp[count] = slope_mean[i]
##            relief_temp[count] = relief_mean[i]
##            Lh_temp[count] = Lh_mean[i]
##            count = count + 1
##        else: print i
##
##    seg = seg_temp[0:count]
##    cht_mean = cht_temp[0:count]
##    slope_mean = slope_temp[0:count]
##    relief_mean = relief_temp[0:count]
##    Lh_mean = Lh_temp[0:count]

    print " done."
    print "\tWriting to " + meanfile + "...",
    f = open(meanfile,'w')
    f.write('seg seg cht_mean slope_mean relief_mean Lh_mean')
    for i in range (0,max_seg):
        f.write('\n' + str(i) + " ")
        f.write(str(seg[i]) + " " + str(cht_mean[i]) + " " + str(slope_mean[i]) + " " + str(relief_mean[i]) + " " + str(Lh_mean[i]))
            

    f.close()

    print " done.\n\n"

if __name__ == "__main__":
    infile = 'db_hillslope_metrics.txt' #filename for raw output of hillslope_tracing.out
    mean_file = 'db_mean_metrics.txt'	#filename for mean metrics for each hilltop segment
    get_means(infile, mean_file)		#launcg
