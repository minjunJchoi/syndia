import imp;
try:
    imp.find_module('MDSplus');
    found = True;
    import MDSplus as MDS;
except ImportError:
    found = False;

import MDSplus as MDS
mds = MDS.Connection("172.17.250.25:8005")
mds.openTree("EFITRT1",18597)
mds.closeTree("EFITRT1",18597)
eq = mds.openTree("EFITRT1",18597)
bcentr=mds.get("\\bcentr")
gtime = mds.get("\\gtime")
def run(shot,time_i,treename='EFITRT1'):
   import os, sys
   import numpy as np
   import matplotlib.pyplot as plt
     
   PARAMETERS=['\\bcentr','\\bdry','\\cpasma','\\epoten','\\ffprim',
               '\\fpol','\\gtime','\\lim','\\limitr','\\mh','\\mw','\\nbdry',
               '\\pprime','\\pres','\\psin','\\psirz','\\qpsi','\\r','\\rgrid1',
               '\\rhovn','\\rmaxis','\\rzero','\\ssibry','\\ssimag','\\xdim','\\z',
               '\\zdim','\\zmaxis','\\zmid'];
   
   geqdsk=[];

   #treename='EFITRT1'
   #treename='EFIT01'
   
   #for efit01 tree

   with MDS.Connection(server="172.17.100.200:8005") as mds:
      try:
          eq=mds.openTree(treename,shot);
      except: 
          print "Error #1"
      else:
          try:
              for signame in PARAMETERS:
                  print 'reading ...',signame
                  temp = mds.get(signame).data();
                  geqdsk.append(temp); 
                  #print geqdsk[PARAMETERS.index(signame)];
              
          except:
              print "Can not reading the signal\n Quit the program";          
              sys.exit(0);
          else:
              print "END of reading"
              #plt.show();


   mds.close(treename,shot);
   
   index_time = PARAMETERS.index('\\gtime');

   if treename == 'EFITRT1':
       pass;
   else:
       geqdsk[index_time]/=1.0e3;
   len_time = len(geqdsk[index_time])
   DTmin = 100.;
   for i in range(len_time):
       DT = time_i - geqdsk[index_time][i]
       if  DT <=DTmin:
           DTmin = DT;
           t_index = i;

   print("Time what we use is %10.5f\n"%geqdsk[index_time][t_index]);
   

   ###### Writing to file #######
   for i in [t_index]:
       time_fileout = time*1000;
       print '%06d'%time_fileout,time
       file_name='kstar_%s_%05d_%06d.geqdsk'%(treename,shot,time_fileout);
       print 'writing..',file_name
       
       f=open(file_name,"w");

       nw=geqdsk[PARAMETERS.index('\\mw')][i];
       nh=geqdsk[PARAMETERS.index('\\mh')][i];

       rleft  = geqdsk[PARAMETERS.index('\\rgrid1')][i]; 
       rdim   = geqdsk[PARAMETERS.index('\\xdim')][i]; 
       rright = rdim+rleft;
       rcentr = geqdsk[PARAMETERS.index('\\rzero')][i]; 
       zdim   = geqdsk[PARAMETERS.index('\\zdim')][i]; 
       zmid   = geqdsk[PARAMETERS.index('\\zmid')][i];

       rmaxis = geqdsk[PARAMETERS.index('\\rmaxis')][i];
       zmaxis = geqdsk[PARAMETERS.index('\\zmaxis')][i];

       simag  = geqdsk[PARAMETERS.index('\\ssimag')][i]; 
       sibry  = geqdsk[PARAMETERS.index('\\ssibry')][i];

       bcentr = geqdsk[PARAMETERS.index('\\bcentr')][i];
       current= geqdsk[PARAMETERS.index('\\cpasma')][i];

       header = "%s  #%08d%08dms                       0"%(treename,shot,time_fileout)+str("%4i"%nw)+str("%4i"%nh)+"\n";
       f.write(header);

       f.write(str("%16.9e"%rdim)+str("%16.9e"%zdim)+str("%16.9e"%rcentr)+str("%16.9e"%rleft)+str("%16.9e"%zmid)+"\n");
       f.write(str("%16.9e"%rmaxis)+str("%16.9e"%zmaxis)+str("%16.9e"%simag)+str("%16.9e"%sibry)+str("%16.9e"%bcentr)+"\n");
       f.write(str("%16.9e"%current)+str("%16.9e"%simag)+str("%16.9e"%0)+str("%16.9e"%rmaxis)+str("%16.9e"%0)+"\n");
       f.write(str("%16.9e"%zmaxis)+str("%16.9e"%0)+str("%16.9e"%sibry)+str("%16.9e"%0)+str("%16.9e"%0)+"\n");

       # profile 

       write_profile(f,geqdsk[PARAMETERS.index('\\fpol')][i]);
       write_profile(f,geqdsk[PARAMETERS.index('\\pres')][i]);
       write_profile(f,geqdsk[PARAMETERS.index('\\ffprim')][i]);
       write_profile(f,geqdsk[PARAMETERS.index('\\pprime')][i]);

       #2D psi ...

       l=0;
       for w in range(nw):
          for h in range(nh):
             f.write(str("%16.9e"%geqdsk[PARAMETERS.index('\\psirz')][i][w][h]));
             l=l+1;
             if(l==5):
                f.write(str("\n"));
                l=0;
            
       if(l>0): f.write(str("\n"))

       # qprofile

       write_profile(f,geqdsk[PARAMETERS.index('\\qpsi')][i]);

       # bdry
       
       nbdry = geqdsk[PARAMETERS.index('\\nbdry')][i];
       nlimt = geqdsk[PARAMETERS.index('\\limitr')][i];
       
       f.write(str("%4i"%nbdry)+str("%4i"%nlimt)+str("\n"));

       l=0         
       for w in range(nbdry):             
         f.write(str("%16.9e"%geqdsk[PARAMETERS.index('\\bdry')][i][w][0]));
         l=l+1;
         if(l==5):
                f.write(str("\n"));
                l=0;
         f.write(str("%16.9e"%geqdsk[PARAMETERS.index('\\bdry')][i][w][1]));
         l=l+1;
         if(l==5):
                f.write(str("\n"));
                l=0;
       if(l>0): f.write(str("\n"))

       l=0         
       for w in range(nlimt):
         f.write(str("%16.9e"%geqdsk[PARAMETERS.index('\\lim')][w][0]));
         l=l+1;
         if(l==5):
                f.write(str("\n"));
                l=0;
         f.write(str("%16.9e"%geqdsk[PARAMETERS.index('\\lim')][w][1]));
         l=l+1;
         if(l==5):
                f.write(str("\n"));
                l=0;
       if(l>0): f.write(str("\n"))

       
       f.close();
        
 
