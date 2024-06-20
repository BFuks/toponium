#!/usr/bin/env python3

# imports
import gzip;
import math;
import os;
import re;
import sys;


#safety 
if not os.path.exists('unweighted_events.lhe.gz'):
    print("    --> Error, cannot find the event file");
    sys.exit();

# IO files
in__file = gzip.open('unweighted_events.lhe.gz', 'r');
out_file = gzip.open('unweighted_events_new.lhe.gz','wb');

# Reading variables initialisation
in_event    = False;
event = []

# colour patterns:
b_pattern = re.compile(br'^\s*5+\s+\d+\s+\d+\s+\d+\s+(\d+)\s')
bbar_pattern = re.compile(br'^\s*-5\s+\d+\s+\d+\s+\d+\s+\d+\s+(\d+)\s')
t_pattern = re.compile(br'^\s*6+\s+\d+\s+\d+\s+\d+\s+(\d+)\s')
tbar_pattern = re.compile(br'^\s*-6\s+\d+\s+\d+\s+\d+\s+\d+\s+(\d+)\s')

# Momenta patterns
pt_pattern = re.compile(br'^\s*6\s+\d+\s+\d+\s+\d+\s+\d+\s+0\s+([+-]?\d+\.\d+e[+-]?\d+)\s+([+-]?\d+\.\d+e[+-]?\d+)\s+([+-]?\d+\.\d+e[+-]?\d+)\s+([+-]?\d+\.\d+e[+-]?\d+)\s')
ptbar_pattern = re.compile(br'^\s*-6\s+\d+\s+\d+\s+\d+\s+0\s+\d+\s+([+-]?\d+\.\d+e[+-]?\d+)\s+([+-]?\d+\.\d+e[+-]?\d+)\s+([+-]?\d+\.\d+e[+-]?\d+)\s+([+-]?\d+\.\d+e[+-]?\d+)\s')

# Reading the file
for line in in__file:
    # Initialisation of a new event
    if b'<event>' in line:
        in_event = True;
        event.append(line);
        continue;

    # getting the number of particles in the event and adding one unit (for the eta state)
    elif in_event:
        event.append(line);

        # Closing the event
        if b'</event>' in line:
           in_event = False;

           # get the antiquark colours
           tbar_colour_code = next((match.group(1).decode() for line in event if (match := tbar_pattern.match(line))), None)
           bbar_colour_code = next((match.group(1).decode() for line in event if (match := bbar_pattern.match(line))), None)

           # get the quark colours
           t_colour_code = next((match.group(1).decode() for line in event if (match := t_pattern.match(line))), None)
           b_colour_code = next((match.group(1).decode() for line in event if (match := b_pattern.match(line))), None)
           if b_colour_code!=t_colour_code or bbar_colour_code!=tbar_colour_code:
               print('ERROR:', event);
               sys.exit()

           # get the toponium momentum
           pt = next((match.groups() for line in event if (match := pt_pattern.match(line))), None)
           ptbar = next((match.groups() for line in event if (match := ptbar_pattern.match(line))), None)
           # Convert bytes to float
           pt = tuple(float(x) for x in pt)
           ptbar = tuple(float(x) for x in ptbar)
           peta = tuple(a + b for a, b in zip(pt, ptbar))

           # updating the event
           ninit=0
           for entry in event:
               # antitop or antibottom
               if entry.decode().strip().startswith("-6") or entry.decode().strip().startswith("-5"):
                   out_file.write(entry.decode().replace(bbar_colour_code,'505').encode())
               # top or bottom
               elif entry.decode().strip().startswith("6") or entry.decode().strip().startswith("5"):
                   out_file.write(entry.decode().replace(b_colour_code,'505').encode())
               # gluon
               elif ' -1 ' in entry.decode():
                   out_file.write(entry.decode().replace(b_colour_code,'506').replace(bbar_colour_code,'506').encode())
                   ninit+=1
                   if ninit==2:
                       toponium = ['      ', 32, '', 2,'  ',1,'  ', 2,'  ', 0,'  ', 0,
                           f"+{peta[0]:.11e}".replace('+-','-'),
                           f"+{peta[1]:.11e}".replace('+-','-'),
                           f"+{peta[2]:.11e}".replace('+-','-'),
                           f"+{peta[3]:.11e}".replace('+-','-'), 
                           f"+{math.sqrt(peta[3]**2-peta[2]**2-peta[1]**2-peta[0]**2):.11e}", f"{0:.4e}", f"{0:.4e}\n"];
                       toponium = [str(x).encode() for x in toponium];
                       out_file.write(" ".encode().join(toponium))
               # other stuff
               else:
                   out_file.write(entry)
           event = [];
        continue

    # outside a given event -> copy the line
    if not in_event:
        out_file.write(line);
        continue

    # else
    print('ERROR = ', line);
    sys.exit();


# closing files
in__file.close();
out_file.close();

# copy files
#os.system('mv unweighted_events.lhe.gz unweighted_events_old.lhe.gz');
#os.system('mv unweighted_events_new.lhe.gz unweighted_events.lhe.gz');

