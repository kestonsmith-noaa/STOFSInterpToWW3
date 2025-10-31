import os
import argparse
import time
import numpy as np
import netCDF4 as nc
import math
import sys
from scipy.interpolate import RegularGridInterpolator
import re
import ComputeMeshToMeshInterpWeights as mshint
 
def loadWW3MeshCoords(fl):
    f=open(fl, 'r')
    header = f.readline() 
    header = f.readline() 
    header = f.readline() 
    header = f.readline() 
    header = f.readline() # number of nodes
    nn = re.findall(r'\d+', header)
    nn=int(nn[0])
    print(nn)
    xi=np.zeros(nn)
    yi=np.zeros(nn)
    k=0
    for i in range(nn):
        A = f.readline()
        #print(A)
        values = A.split(" ")
#        print(values)
        xi[k]=float(values[2])
        yi[k]=float(values[4])
        k=k+1
    return xi, yi 
    f.close


def WriteInterpJobscript(fl,N, ComputeNodes):
    with open(fl, 'w') as f:
        #yi[k]=float(values[4])
        f.write("#!/bin/bash \n")
        f.write("#SBATCH --job-name=STOFS_interp_masterscript \n")
#        f.write("#SBATCH --ntasks="+str(N)+" \n")
        f.write("#SBATCH --ntasks=1 \n") # ntasks per interpolation
        f.write("#SBATCH --time=08:00:00 \n") 
        f.write("#SBATCH --output=mpi_test_%j.log \n")
        f.write("#SBATCH --error=%j.err \n")
        f.write("#SBATCH --account=marine-cpu \n")
        f.write("#SBATCH --nodes="+str(ComputeNodes)+" \n")
        f.write("#SBATCH --ntasks-per-core=1"+" \n")
        f.write("#SBATCH --array=0-"+str(N-1)+" \n")

        f.write(" \n")

        f.write("module purge \n")
        f.write("module use /scratch4/NCEPDEV/marine/Ali.Salimi/Hera_Data/HR4-OPT/FromJessica/Keston/ICunstructuredRuns15km-implicit-450s/global-workflow/sorc/ufs_model.fd/modulefiles \n")
        f.write("module load ufs_ursa.intel \n")
        f.write("module load py-scipy/1.14.1 \n")
        f.write("module load py-netcdf4/1.7.1.post2 \n")
        f.write("pip list \n")
        f.write("srun python3  InterpSTOFS2mesh.part.py $SLURM_ARRAY_TASK_ID > InterpJob.$SLURM_ARRAY_TASK_ID.out")
#        for k in range(N):
#            f.write("srun -n 1 python3 InterpSTOFS2mesh.part.py "+str(k)+" > InterpJob."+str(k)+".out & \n")
        f.write("wait\n")
        f.write("##run this after the full job array is compleate, need to test")
        f.write("python3 ConvertOutput2Netcdf.py 50")





mesh="global_1deg_unstr"
OutDir=mesh+".files/"
# Create the output directory------------------------------------------
try:
    os.mkdir(OutDir)
    print(f"Directory '{OutDir}' created successfully.")
except FileExistsError:
    print(f"Directory '{OutDir}' already exists. Proceeding ...")
except PermissionError:
    print(f"Permission denied: Unable to create '{OutDir}'.")

#fl="meshes/global_1deg_unstr.msh"
fl="meshes/"+mesh+".msh"
f=open(fl, 'r')
for k in range(5):
    header = f.readline() 
#nn = re.findall(r'\d+', header)
#nn=int(nn[0])
#print("nn = "+str(nn))
xi, yi =loadWW3MeshCoords(fl)
Nparts=int(sys.argv[1])
nn=xi.shape[0]
NodesPerProc=math.ceil(nn/Nparts)
print(str(NodesPerProc)+ " " + str(nn))
a = range(nn)
NodeList = [] 
for i in range(0, len(a), NodesPerProc):  # Slice list in steps of n
    NodeList.append(a[i:i + NodesPerProc])

print(NodeList)
for k in range(Nparts):
    flo=OutDir+'NodeList.'+str(k)+'.txt'
    f=open(flo, 'w')
    f.write("List of nodes for job "+str(k) + "\n")
    List=NodeList[k]
    List
    n=len(List)
    f.write(str(n) + "\n")
    for j in range(n):
        f.write(str(List[j]) + '\n')
    f.close

#write jobcard to do interpolation
WriteInterpJobscript("jobcardInterp2Mesh",Nparts,1)
