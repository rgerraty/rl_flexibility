#!/bin/bash
#   Copyright (C) 2012 University of Oxford
#
#   Part of FSL - FMRIB's Software Library
#   http://www.fmrib.ox.ac.uk/fsl
#   fsl@fmrib.ox.ac.uk
#   
#   Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
#   Imaging of the Brain), Department of Clinical Neurology, Oxford
#   University, Oxford, UK
#   
#   
#   LICENCE
#   
#   FMRIB Software Library, Release 5.0 (c) 2012, The University of
#   Oxford (the "Software")
#   
#   The Software remains the property of the University of Oxford ("the
#   University").
#   
#   The Software is distributed "AS IS" under this Licence solely for
#   non-commercial use in the hope that it will be useful, but in order
#   that the University as a charitable foundation protects its assets for
#   the benefit of its educational and research purposes, the University
#   makes clear that no condition is made or to be implied, nor is any
#   warranty given or to be implied, as to the accuracy of the Software,
#   or that it will be suitable for any particular purpose or for use
#   under any specific conditions. Furthermore, the University disclaims
#   all responsibility for the use which is made of the Software. It
#   further disclaims any liability for the outcomes arising from using
#   the Software.
#   
#   The Licensee agrees to indemnify the University and hold the
#   University harmless from and against any and all claims, damages and
#   liabilities asserted by third parties (including claims for
#   negligence) which arise directly or indirectly from the use of the
#   Software or the sale of any products based on the Software.
#   
#   No part of the Software may be reproduced, modified, transmitted or
#   transferred in any form or by any means, electronic or mechanical,
#   without the express permission of the University. The permission of
#   the University is not required if the said reproduction, modification,
#   transmission or transference is done without financial return, the
#   conditions of this Licence are imposed upon the receiver of the
#   product, and all original and amended source code is included in any
#   transmitted product. You may be held legally responsible for any
#   copyright infringement that is caused or encouraged by your failure to
#   abide by these terms and conditions.
#   
#   You are not permitted under this Licence to use this Software
#   commercially. Use for which any financial return is received shall be
#   defined as commercial use, and includes (1) integration of all or part
#   of the source code or the Software into a product for sale or license
#   by or on behalf of Licensee to third parties or (2) use of the
#   Software or any derivative of it for research with the final aim of
#   developing software products for sale or license to a third party or
#   (3) use of the Software or any derivative of it for research with the
#   final aim of developing non-software products for sale or license to a
#   third party, or (4) use of the Software to provide any service to an
#   external organisation for which payment is received. If you are
#   interested in using the Software commercially, please contact Isis
#   Innovation Limited ("Isis"), the technology transfer company of the
#   University, to negotiate a licence. Contact details are:
#   innovation@isis.ox.ac.uk quoting reference DE/9564.

##10/03/13-- edited by Raphael Gerraty to take in a 9 parameter file including 
##CSF, white matter, and whole brain signals for use in resting-state confound 
##regression.  

if [ $# -lt 2 ] ; then
    cat << EOF
Usage: mp_diff regparam.dat diffregparam.dat

Custom version of FSL confound regressor script for use with resting-state data.
Creates file with 27 columns; the first 9 are the square of the motion and csf, wb, and wm
parameters, the next 9 are the temporal difference of these parameters, 
and the next 9 are the square of the differenced values.  This is useful for
accounting for 'spin history' effects, and variation not 
otherwise accounted for by confound correction.

\$Id: mp_diffpow.sh,v 1.2 2012/08/30 13:50:50 mwebster Exp $
EOF
    exit 1
fi

f=`echo $2 | sed 's/\....$//'`

cat <<EOF > /tmp/$$-mp-diffpow
{
  if (NR==1) {
    mp1=\$1;mp2=\$2;mp3=\$3;mp4=\$4;mp5=\$5;mp6=\$6;mp7=\$7;mp8=\$8;mp9=\$9;
  }
  printf("  %+0.7e  %+0.7e  %+0.7e  %+0.7e  %+0.7e  %+0.7e  %+0.7e  %+0.7e  %+0.7e    %+0.7e  %+0.7e  %+0.7e  %+0.7e  %+0.7e  %+0.7e  %+0.7e  %+0.7e  %+0.7e    %+0.7e  %+0.7e  %+0.7e  %+0.7e  %+0.7e  %+0.7e  %+0.7e  %+0.7e  %+0.7e\n",
         \$1^2,\$2^2,\$3^2,\$4^2,\$5^2,\$6^2,\$7^2,\$8^2,\$9^2,
         \$1-mp1,\$2-mp2,\$3-mp3,\$4-mp4,\$5-mp5,\$6-mp6,\$7-mp7,\$8-mp8,\$9-mp9,
         (\$1-mp1)^2,(\$2-mp2)^2,(\$3-mp3)^2,(\$4-mp4)^2,(\$5-mp5)^2,(\$6-mp6)^2,(\$7-mp7)^2,(\$8-mp8)^2,(\$9-mp9)^2);
  mp1=\$1;mp2=\$2;mp3=\$3;mp4=\$4;mp5=\$5;mp6=\$6;mp7=\$7;mp8=\$8;mp9=\$9;
}
EOF
awk -f /tmp/$$-mp-diffpow "$1" > ${f}.dat

/bin/rm /tmp/$$-mp-diffpow
