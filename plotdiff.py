#!/Users/clidman/Science/Programs/Ureka/variants/common/bin/python
# Python program to plot difference between poschecks

# 2018/09/27
# Reads in the values using sds routines
# 2018/12/01
# Option to exclude the offset
# 2020/02/14
# Prints out the offset
# 2022/03/28
# MR adapting to run on aatlxy

# To run code on AATLXY, go to /instsoft/2dF/config/poscheck_diff/
# Then run below command,
# python plotdiff.py -a -1 /path/to/old/sdsraw.sds -2 /path/to/new/sdsraw.sds
# PDF document will be saved in /instsoft/2dF/config/poscheck_diff/

import numpy
import matplotlib.pyplot as plt
import subprocess as sp
from optparse import OptionParser
import copy
import warnings

warnings.simplefilter("ignore", FutureWarning)

class nf(float):
    def __repr__(self):
        str = '%.1f' % (self.__float__(),)
        if str[-1] == '0':
            return '%.0f' % self.__float__()
        else:
            return '%.1f' % self.__float__()

class poscheck:
    """A poscheck result"""
    def __inti__(self):
        self.c0=0.0
        self.c1=0.0
        self.c2=0.0
        self.c3=0.0
        self.c4=0.0
        self.c5=0.0
        self.x_off=0.0        # microns
        self.y_off=0.0        # microns
        self.a=13366713.7     # Number of microns per radian
        self.b=1089450000.0
        self.c=621000000000.0
        self.d=331200000000000.0

        return

    def convert(self,x,y):
        # Convert from microns to radians
        x=x/self.a
        y=y/self.a

        # Distortion correction first
        x1=(x-self.x_off/self.a)
        y1=(y-self.y_off/self.a)

        # self.x_off and self.y_off in microns
        # x1 and y1 are in radians
        # x and y are in radians

        r2=x1**2.+y1**2.

        # Interesting notation
        x2=x*(self.a+r2*(self.b+r2*(self.c+r2*self.d)))
        y2=y*(self.a+r2*(self.b+r2*(self.c+r2*self.d)))

        # x2 and y2 are in microns

        # Linear correction second

        x_3=self.c0+self.c1*x2+self.c2*y2
        y_3=self.c3+self.c4*x2+self.c5*y2

        # x_3 and y_3 are in microns

        # Convert from radians to arc seconds
        factor = self.a / (3600 * 180 / numpy.pi)

        return x_3/factor,y_3/factor

parser = OptionParser()

parser.add_option("-1", "--file1", dest="file1",
                  help="The first poscheck file")

parser.add_option("-2", "--file2", dest="file2",
                  help="The second poscheck file")

parser.add_option("-o", "--offset", dest="offset", default=0.0, type=float,
                  help="Offset in x (microns)")

parser.add_option("-s", "--scale", dest="scale", default=0.00005, type=float,
                  help="Offset in x")

parser.add_option("-l", "--length", dest="length", default=1.0, type=float,
                  help="Length of arrow in arc seconds")

parser.add_option("-a", "--applyOffset", dest="applyOffset", default=False,
                  action="store_true",
                  help="Apply offset")


(options, args) = parser.parse_args()

# Add path for sdsdump. Shells do not have drama commands loaded.
# Done in a hack method for ease as this script should not be used elsewhere.
sdspath = '/instsoft/drama/release/sds/r1_4_3_92/linux_x86/'

p1=poscheck()
p2=poscheck()

cmd="sdsdump %s PosCheckResults.newDist" % (options.file1)
cmd = sdspath + cmd
p = sp.Popen(cmd,stdout=sp.PIPE,shell=True)
output = p.communicate()[0].split()

cmd="sdsdump %s PosCheckResults.newCoeffs" % (options.file1)
cmd = sdspath + cmd
p = sp.Popen(cmd,stdout=sp.PIPE,shell=True)
output2 = p.communicate()[0].split()

p1.c0=float(output2[0])
p1.c1=float(output2[1])
p1.c2=float(output2[2])
p1.c3=float(output2[3])
p1.c4=float(output2[4])
p1.c5=float(output2[5])
p1.x_off=float(output[4])
p1.y_off=float(output[5])
p1.a=13366713.699999999255  #1e6
p1.b=1089450000.0           #1e12
p1.c=621000000000.0         #1e18
p1.d=331200000000000.0      #1e24

cmd="sdsdump %s PosCheckResults.newDist" % (options.file2)
cmd = sdspath + cmd
p = sp.Popen(cmd,stdout=sp.PIPE,shell=True)
output = p.communicate()[0].split()

cmd="sdsdump %s PosCheckResults.newCoeffs" % (options.file2)
cmd = sdspath + cmd
p = sp.Popen(cmd,stdout=sp.PIPE,shell=True)
output2 = p.communicate()[0].split()

p2.c0=float(output2[0])
p2.c1=float(output2[1])
p2.c2=float(output2[2])
p2.c3=float(output2[3])
p2.c4=float(output2[4])
p2.c5=float(output2[5])
p2.x_off=float(output[4])
p2.y_off=float(output[5])
p2.a=13366713.699999999255
p2.b=1089450000.0
p2.c=621000000000.0
p2.d=331200000000000.0

if options.offset != 0.0:
    p2=copy.copy(p1)
    p2.x_off+=options.offset

# x and y are in microns
x=numpy.arange(-250000,250000,25000.)
y=numpy.arange(-250000,250000,25000.)

xx,yy=numpy.meshgrid(x,y)

# x_3 and y_3 are in arc seconds
x_3_p1,y_3_p1=p1.convert(xx,yy)
x_3_p2,y_3_p2=p2.convert(xx,yy)

x_offset=x_3_p1-x_3_p2
y_offset=y_3_p1-y_3_p2

length=244957.8

offset=numpy.sqrt((x_offset)**2.+(y_offset)**2.)

# Compute the average over the middle
selection =(xx*xx+yy*yy < (length  * 0.6667)**2.)

fig=plt.figure(figsize=(10,10))
ax=fig.add_subplot(111)



if options.applyOffset:
    levels = numpy.arange(-2, 2.1, 0.1)
    CS=ax.contour(xx,yy,offset-numpy.average(offset[selection]),levels)
else:
    levels = numpy.arange(-2, 2.1, 0.1) + numpy.average(offset[selection])
    CS=ax.contour(xx,yy,offset,levels)

avg_offset = 'Average Offset %4.1f %4.1f (arc seconds)' % (numpy.average(x_offset[selection]),numpy.average(y_offset[selection]))

print(avg_offset)

if options.applyOffset:
    q = ax.quiver(xx, yy, x_offset-numpy.average(x_offset[selection]), y_offset-numpy.average(y_offset[selection]),scale_units='xy',scale=options.scale)
else:
    q = ax.quiver(xx, yy, x_offset, y_offset,scale_units='xy',scale=options.scale)

scaling=0.005 * options.length # the 0.005 factor is used to get the scale roughly correct. It is unclear where this comes from

#ax.quiverkey(q, X=0.9, Y=1.05, U=scaling,
#             label='%4.1f" Offset' % options.length, labelpos='E')

#q = ax.quiver(xx, yy, x_offset, y_offset,scale_units='xy',scale=0.00005)
#ax.quiverkey(q, X=0.3, Y=1.05, U=2,
#             label='1" Offset', labelpos='E')


# Label levels with specially formatted floats
if plt.rcParams["text.usetex"]:
    fmt = r'%r'
else:
    fmt = '%r'

# Recast levels to new class
CS.levels = [nf(val) for val in CS.levels]

ax.clabel(CS, CS.levels, inline=True, fmt=fmt, fontsize=30)


circle1 = plt.Circle((0, 0), length, color='k', fill=False)

ax.add_artist(circle1)
ax.set_xlim(-1*length,length)
ax.set_ylim(-1*length,length)

ax.set_title(avg_offset)

plt.savefig("%s_%s.pdf" % (options.file1.replace(".sds","").replace("../","").replace("/","_"),options.file2.replace(".sds","").replace("../","").replace("/","_")))
plt.show()
plt.close()
