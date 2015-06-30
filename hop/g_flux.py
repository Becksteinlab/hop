import MDAnalysis
import numpy as np
import MDAnalysis.KDTree.NeighborSearch as NS
import copy
import time

def setup(topology,trajectory,water_label,protein_label):

    u=MDAnalysis.Universe(topology,trajectory)
    w=u.selectAtoms('name ' + water_label)
    box=u.atoms.bbox()
    l=u.selectAtoms('resname ' + protein_label)
    p=u.selectAtoms('protein')
    protein_box=p.bbox()
    return [u,w]

def track_counts(topology,trajectory,water_label,protein_label,in_top,in_bottom,write_out_time=1000,timestep=0.001,use_cutoff=False,time_cutoff=1000,filename='rates'):

    universe,water=setup(topology,trajectory,
    water_label,protein_label)
    u=universe
    total_time=universe.trajectory.totaltime
    tracking_list=np.zeros(len(water))
    tracking_list_last=np.zeros(len(water))
    counts_up=0
    counts_down=0
    rate_observation=0
    time_elapsed=0
    crossed_particles_up=[]
    crossed_particles_down=[]
    rates=open(filename+'.txt','w')
    rates.write('timestep' + ' ' + str(timestep) + ' ' + 'upper boundary' + ' ' + str(in_top) + ' ' + 'lower boundary' + ' ' + str(in_bottom) + '\n')
    cumulative_counts=0
    total_residence_time=0
    weighted_residence_time=0
    if use_cutoff==False:
        time_cutoff=total_time
    for ts in universe.trajectory:
            inst_counts_down=0
            inst_counts_up=0
            time_elapsed+=1
            water_z=water.coordinates()[:,2]
            progress=np.divide(universe.trajectory.frame,universe.trajectory.numframes,dtype=np.float64)
            tracking_list_last[:]=tracking_list
            if time_elapsed>=write_out_time:
                time_elapsed=0
            for water_index in xrange(len(water_z)):
                if water_z[water_index]>=in_top:
                    tracking_list[water_index]=1
                    if tracking_list_last[water_index]==-2:
                        counts_up+=1
                elif water_z[water_index]<in_top and water_z[water_index]>in_bottom:
                    if tracking_list_last[water_index]==1 or tracking_list_last[water_index]==2:
                        tracking_list[water_index]=2
                    elif tracking_list_last[water_index]==-1 or tracking_list_last[water_index]==-2:
                        tracking_list[water_index]=-2
                elif water_z[water_index]<in_bottom:
                    tracking_list[water_index]=-1
                    if tracking_list_last[water_index]==2:
                        counts_down+=1
            if universe.trajectory.frame!=0:
                rate_up=np.divide(counts_up,u.trajectory.frame*timestep,dtype=np.float64)
                rates.write(str(rate_up) + '\n')
                rate_down=np.divide(counts_down,u.trajectory.frame*timestep,dtype=np.float64)
                rates.write(str(rate_down) + '\n')
                
                rate=np.divide(counts_up-counts_down,u.trajectory.frame*timestep,dtype=np.float64)
                inst_rate=np.divide(inst_counts_up-inst_counts_down,u.trajectory.frame*timestep,dtype=np.float64)
    return (rate_up,rate_down,crossed_particles_up,crossed_particles_down)

