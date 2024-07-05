#!/usr/bin/env python3

# imports
import gzip;
import math;
import os;
import re;
import sys;

# Debug 
debug = False

#safety 
if not os.path.exists('unweighted_events.lhe.gz'):
    print("    --> Error, cannot find the event file");
    sys.exit();

# IO files
in__file = gzip.open('unweighted_events.lhe.gz', 'r');
out_file = gzip.open('unweighted_events_new.lhe.gz','wb');

# Reading variables initialisation
in_event     = False;
missing_tops = 0
event        = []

# colour patterns:
b_pattern = re.compile(br'^\s*5+\s+\d+\s+\d+\s+\d+\s+(\d+)\s')
bbar_pattern = re.compile(br'^\s*-5\s+\d+\s+\d+\s+\d+\s+\d+\s+(\d+)\s')
t_pattern = re.compile(br'^\s*6+\s+\d+\s+\d+\s+\d+\s+(\d+)\s')
tbar_pattern = re.compile(br'^\s*-6\s+\d+\s+\d+\s+\d+\s+\d+\s+(\d+)\s')

# Momenta patterns
pt_pattern = re.compile(br'^\s*6\s+\d+\s+\d+\s+\d+\s+\d+\s+0\s+([+-]?\d+\.\d+e[+-]?\d+)\s+([+-]?\d+\.\d+e[+-]?\d+)\s+([+-]?\d+\.\d+e[+-]?\d+)\s+([+-]?\d+\.\d+e[+-]?\d+)\s')
ptbar_pattern = re.compile(br'^\s*-6\s+\d+\s+\d+\s+\d+\s+0\s+\d+\s+([+-]?\d+\.\d+e[+-]?\d+)\s+([+-]?\d+\.\d+e[+-]?\d+)\s+([+-]?\d+\.\d+e[+-]?\d+)\s+([+-]?\d+\.\d+e[+-]?\d+)\s')

def updatemothers(entry, pid=''):
    myentry = entry.decode().strip().split()
    if pid=='':
        myentry[2] = str(int(myentry[2])+1)
        myentry[3] = str(int(myentry[3])+1)
    else:
        myentry[2] = str(pid)
        myentry[3] = str(pid)
    myentry  = "       ".encode() + (" ".join(myentry)+"\n").encode()
    return myentry

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

           # Do we have a missing top ?
           pdgcodes = [int(re.match(r'\s*(-?\d+)', entry.decode()).group()) for entry in event if re.match(r'\s*(-?\d+)', entry.decode())]
           if not 6 in pdgcodes:
              missing_tops+=1
              pos_w, prt_w = next(((i, line.decode().strip()) for i, line in enumerate(event) if line.decode().strip().startswith('24')), (None, None))
              prt_w = "       "+prt_w.replace("1    2", str(pos_w-1)+"    "+str(pos_w-1))
              prt_b = "       "+next((line.decode().strip() for line in event if line.decode().strip().startswith('5')), None).replace("1    2", str(pos_w-1)+"    "+str(pos_w-1))
              b_colour_code = next((match.group(1).decode() for line in event if (match := b_pattern.match(line))), None)
              event[pos_w + 1:] = [updatemothers(line) if '>' not in line.decode() and int(line.decode().split()[2]) >= (pos_w - 1) else line for line in event[pos_w + 1:]]
              event = [prt_w.encode() + b'\n' if line.decode().strip().startswith('24') else line for line in event]
              event = [prt_b.encode() + b'\n' if line.decode().strip().startswith('5') else line for line in event]
              prt_t = [str(float(prt_w.split()[i]) + float(prt_b.split()[i])) for i in range(6,10)]
              prt_t = " ".join(["        6  2    1    2  "+ str(b_colour_code)+"    0 "] + \
                  [f"+{float(prt_t[0]):.11e}".replace('+-','-')] + \
                  [f"+{float(prt_t[1]):.11e}".replace('+-','-')] + \
                  [f"+{float(prt_t[2]):.11e}".replace('+-','-')] + \
                  [f"+{float(prt_t[3]):.11e}".replace('+-','-')] + \
                  [f"+{math.sqrt(float(prt_t[3])**2 - float(prt_t[1])**2 - float(prt_t[2])**2 - float(prt_t[0])**2):.11e}"] + \
                  ["0.0000e+00 0.0000e+00\n"])
              event.insert(pos_w, prt_t.encode())

           if not -6 in pdgcodes:
              missing_tops+=1
              pos_w, prt_w = next(((i, line.decode().strip()) for i, line in enumerate(event) if line.decode().strip().startswith('-24')), (None, None))
              prt_w = "       "+prt_w.replace("1    2", str(pos_w-1)+"    "+str(pos_w-1))
              prt_b = "       "+next((line.decode().strip() for line in event if line.decode().strip().startswith('-5')), None).replace("1    2", str(pos_w-1)+"    "+str(pos_w-1))
              b_colour_code = next((match.group(1).decode() for line in event if (match := bbar_pattern.match(line))), None)
              event[pos_w + 1:] = [updatemothers(line) if '>' not in line.decode() and int(line.decode().split()[2]) >= (pos_w - 1) else line for line in event[pos_w + 1:]]
              event = [prt_w.encode() + b'\n' if line.decode().strip().startswith('-24') else line for line in event]
              event = [prt_b.encode() + b'\n' if line.decode().strip().startswith('-5') else line for line in event]
              prt_t = [str(float(prt_w.split()[i]) + float(prt_b.split()[i])) for i in range(6,10)]
              prt_t = " ".join(["       -6  2    1    2  0    "+str(b_colour_code)+" "] + \
                  [f"+{float(prt_t[0]):.11e}".replace('+-','-')] + \
                  [f"+{float(prt_t[1]):.11e}".replace('+-','-')] + \
                  [f"+{float(prt_t[2]):.11e}".replace('+-','-')] + \
                  [f"+{float(prt_t[3]):.11e}".replace('+-','-')] + \
                  [f"+{math.sqrt(float(prt_t[3])**2 - float(prt_t[1])**2 - float(prt_t[2])**2 - float(prt_t[0])**2):.11e}"] + \
                  ["0.0000e+00 0.0000e+00\n"])
              event.insert(pos_w, prt_t.encode())

           # get the antiquark colours
           tbar_colour_code = next((match.group(1).decode() for line in event if (match := tbar_pattern.match(line))), None)
           bbar_colour_code = next((match.group(1).decode() for line in event if (match := bbar_pattern.match(line))), None)

           # get the quark colours
           t_colour_code = next((match.group(1).decode() for line in event if (match := t_pattern.match(line))), None)
           b_colour_code = next((match.group(1).decode() for line in event if (match := b_pattern.match(line))), None)
           if b_colour_code!=t_colour_code or bbar_colour_code!=tbar_colour_code:
               print('ERROR with colour codes:')
               print('     -> tbar = ,', tbar_colour_code)
               print('     -> bbar = ,', bbar_colour_code)
               print('     -> t    = ,', t_colour_code)
               print('     -> b    = ,', b_colour_code)
               print(" ".join([x.decode() for x in event if not '<' in x.decode()]))
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
           init_line_written=False
           for entry in event:
               # init line (increment of the number of particle)
               if not init_line_written and not entry.decode().strip().startswith("<event>"):
                   init_line_written = True
                   myentry = entry.decode().strip().split()
                   myentry[0]=str(int(myentry[0])+1+missing_tops)
                   out_file.write((" ".join(myentry)+'\n').encode())
                   if debug: print(" ".join(myentry)+'\n')
               # antitop or antibottom
               elif entry.decode().strip().startswith("-6"):
                   out_file.write(updatemothers(entry, pid="3").decode().replace(bbar_colour_code,'505').encode())
                   if debug: print(updatemothers(entry, pid="3").decode().replace(bbar_colour_code,'505'))
               elif entry.decode().strip().startswith("-5"):
                   out_file.write(updatemothers(entry).decode().replace(bbar_colour_code,'505').encode())
                   if debug: print(updatemothers(entry).decode().replace(bbar_colour_code,'505'))
               # top or bottom
               elif entry.decode().strip().startswith("6"):
                   out_file.write(updatemothers(entry, pid="3").decode().replace(b_colour_code,'505').encode())
                   if debug: print(updatemothers(entry, pid="3").decode().replace(b_colour_code,'505'))
               elif entry.decode().strip().startswith("5"):
                   out_file.write(updatemothers(entry).decode().replace(b_colour_code,'505').encode())
                   if debug: print(updatemothers(entry).decode().replace(b_colour_code,'505'))
               # gluon
               elif ' -1 ' in entry.decode():
                   out_file.write(entry.decode().replace(b_colour_code,'506').replace(bbar_colour_code,'506').encode())
                   if debug: print(entry.decode().replace(b_colour_code,'506').replace(bbar_colour_code,'506'))
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
                       if debug: print((" ".encode().join(toponium)).decode())
               # other stuff
               else:
                   if any(word in entry.decode() for word in ['wgt', 'rwt', 'rscale', 'totfact', 'event']):
                       out_file.write(entry)
                   else:
                       out_file.write(updatemothers(entry))
                       if debug: print(updatemothers(entry).decode())
           event = [];
           missing_tops = 0
           init_line_written=False
           if debug:
               sys.exit()
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
os.system('mv unweighted_events.lhe.gz unweighted_events_old.lhe.gz');
os.system('mv unweighted_events_new.lhe.gz unweighted_events.lhe.gz');

