from pymol import cmd,stored
from math import sqrt
import numpy as np

class NMRpymol:

    def __init__(self, datafile, sele):
        self.data = np.loadtxt(datafile)
        self.res  = self.data[:,0]
        self.vals = self.data[:,1]
     
        self.myspace = {'bfactorlist': [] }
        cmd.iterate('%s  ' % ( sele ), 'bfactorlist.append((resi, name, b, color))', space=self.myspace )    
        self.sele_resis=[int(i[0]) for i in self.myspace['bfactorlist'] ]


def alterBfactors(  sele, filename, mini=0, maxi=1):
    """alters the bfactors for the specified selection based on a datafile
        within the selected range. Only values for residues within the pymol
        selection are considered for normalisation. Default mini=0, maxi=1.0
    
        The datafile hs the format: 
        residue value

        example usage

        run path/nmrpymol.py

        alterBfactors all, path_to_filename
        alterfactors resid 40:50, path_to_filename, mini=0, maxi=1

    """
    nmrpymol=NMRpymol(filename, sele)
    
    cmd.alter('all' , 'b="0.000"')
    mini, maxi =float(mini), float(maxi)

    #select only the datavalues of the pymol selected residues 
    sele_res, sele_vals= [], []
    for i, r in enumerate(nmrpymol.res):
        if r in nmrpymol.sele_resis:
            sele_res.append(r)
            sele_vals.append(nmrpymol.vals[i])

    xmin, xmax = np.min(sele_vals), np.max(sele_vals)

    data_dict={}    
    for i, r in enumerate(sele_res): 
        if r in nmrpymol.sele_resis:
            norm = (sele_vals[i] - xmin)/(xmax - xmin)
            data_dict['%i'%r]= "%.3f" %(mini + (maxi-mini)*norm)


    for i in nmrpymol.myspace['bfactorlist']:
        try:
            newBfact=float(data_dict[str(int(i[0]))])
        except:
            newBfact=0.000
        
        cmd.alter('resi %i' %(int(i[0])), 'b="%.3f"'%newBfact)
        print('resi %i new b_fact: %.3f' %(int(i[0]),  newBfact ) )
    
    
def color_by_data(  sele, filename, offset=0, colormap='gray70 yellow orange red'):
    """colors a selection based on a datafile using the specified colors in colormap
        supports only 3 [<1, 1<v<2, >2 std] or 4 colors [<0.5, 0.5<v<1, 1<v<2, >2 std]. 
        A residue offset between the datafile and the resi
        
        The datafile hs the format: 
        residue value

        example usage

        run path/nmrpymol.py

        color_by_data 6u3r,  data.txt ,  offset=0, colormap=gray70 yellow red
        color_by_data 6u3r and resid 1:60,  data.txt ,   colormap=gray70 yellow orange red
    """
    nmrpymol=NMRpymol(filename, sele)
    print(nmrpymol.sele_resis)

    colors=colormap.split(' ')
    cmd.color('gray70', sele, quiet=1)
    
    miss_sele=''
    for pdb_res in nmrpymol.sele_resis:
        if pdb_res not in nmrpymol.res:
            print('Missing residue: %i' %pdb_res)
            miss_sele += 'resid %i or ' %pdb_res
    cmd.select('Missing_resis', miss_sele[:-4])
    cmd.disable('Missing_resis')
    
    if len(colors)>4 or len(colors)<3:
        print('Error: Currently only 3 or 4 colors are supported')
        return

    for resid, val in zip(nmrpymol.res, nmrpymol.vals):
        resid=int(resid)+int(offset)
    
        if resid in nmrpymol.sele_resis:
            if len(colors)==4:
                lims=[0.5, 1, 2]
                if val < lims[0]*np.std(nmrpymol.vals):
                    cmd.color(colors[0], '%s and resid %i' %(sele, resid))
                    print('resid %i less than 0.5*stdev (%.3f) color: %s'%(resid, val, colors[0]))
            
                elif val > lims[0]*np.std(nmrpymol.vals) and val < lims[1]*np.std(nmrpymol.vals):
                    cmd.color(colors[1], '%s and resid %i' %(sele, resid))
                    print('resid %i more than 0.5*stdev (%.3f) color: %s'%(resid, val, colors[1]))
            
                elif val > lims[1]*np.std(nmrpymol.vals) and val < lims[2]*np.std(nmrpymol.vals):
                    cmd.color(colors[2], '%s and resid %i' %(sele, resid)) 
                    print('resid %i more than 1.0*stdev (%.3f) color: %s'%(resid, val, colors[2]))
                
                elif val > lims[-1]*np.std(nmrpymol.vals):
                    cmd.color(colors[3], '%s and resid %i' %(sele, resid))
                    print('resid %i more than 2.0*stdev (%.3f) color: %s'%(resid, val, colors[-1]))
                
        
            if len(colors)==3:  
                lims=[1, 2]         
                if val < lims[0]*np.std(nmrpymol.vals):
                    cmd.color(colors[0], '%s and resid %i' %(sele, resid))
                    print('resid %i less than 1.0*stdev (%.3f) color: %s'%(resid, val, colors[0]))
            
                elif val > lims[0]*np.std(nmrpymol.vals) and val < lims[1]*np.std(nmrpymol.vals):
                    cmd.color(colors[1], '%s and resid %i' %(sele, resid))
                    print('resid %i more than 1.0*stdev (%.3f) color: %s'%(resid, val, colors[1]))
                
                elif val > lims[-1]*np.std(nmrpymol.vals):
                    cmd.color(colors[2], '%s and resid %i' %(sele, resid))    
                    print('resid %i more than 2.0*stdev (%.3f) color: %s'%(resid, val, colors[-1]))
                            
        else:
            cmd.color('gray30', '%s and resid %i' %(sele, resid))
    
cmd.extend("alterBfactors", alterBfactors)
cmd.extend("color_by_data", color_by_data)
