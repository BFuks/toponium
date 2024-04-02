#!/usr/bin/env python3

# imports
import gzip;
import math;
import os;
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
n_particles = -1;
prt_counter =  0;

# Reading the file
for line in in__file:
    # Initialisation of a new event
    if b'<event>' in line:
        in_event = True;
        particles = [];
        out_file.write(line);
        continue;

    # getting the number of particles in the event and adding one unit (for the eta state)
    if in_event and n_particles == -1:
        myline = line.split();
        n_particles = int(myline[0]);
        myline[0]=str(n_particles+1).encode();
        out_file.write(' '.encode().join(myline) + '\n'.encode());
        prt_counter = 0;
        continue;

    # we have a particle
    if n_particles > prt_counter:
        prt_counter+=1;
        myline = line.split();
        if   prt_counter==1 and int(myline[0])==21:
            myline[4]  = '501'.encode();
            myline[5]  = '502'.encode();
        elif prt_counter==2 and int(myline[0])==21:
            myline[4]  = '502'.encode();
            myline[5]  = '501'.encode();
        elif int(myline[0]) in [6,5]:
            myline[4]  = '503'.encode();
        elif int(myline[0]) in [-6,-5]:
            myline[5]  = '503'.encode();
        particles.append(myline);
        continue;

    # End of the event
    if n_particles == prt_counter:
        # computing the momentum of the toponium state and updating the decay chain
        px=0.; py=0.; pz=0.; E=0.;
        for index in range(len(particles)):
            pdg = int(particles[index][0]);
            if pdg in [5,-5,24,-24]:
                px += float(particles[index][6]);
                py += float(particles[index][7]);
                pz += float(particles[index][8]);
                E  += float(particles[index][9]);
            if index>1:
                if int(particles[index][2])==1 and int(particles[index][3])==2:
                    particles[index][2]='3'.encode();
                    particles[index][3]='3'.encode();
                elif particles[index][2]==particles[index][3]:
                    particles[index][2] = str(int(particles[index][2])+1).encode();
                    particles[index][3] = str(int(particles[index][3])+1).encode();
        toponium = [32, 2, 1, 2, 0, 0, px, py, pz, E, math.sqrt(E**2-px**2-py**2-pz**2), 0.0000e+00, 0.00000e+0];
        toponium = [str(x).encode() for x in toponium];
        particles.insert(2,toponium);
        for part in particles: out_file.write(' '.encode().join(part) + '\n'.encode());

        # reseting all variables
        in_event = False;
        n_particles = -1;
        prt_counter = 0;
        particles   = [];
        out_file.write(line);
        continue;

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

