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
    #tracking_list=np.zeros((len(water),len(water)))
    tracking_list=np.zeros(len(water))
    #tracking_list_last=np.zeros(len(water))
    tracking_list_last=np.zeros(len(water))
    counts_up=0
    counts_down=0
    rate_observation=0
    time_elapsed=0
    crossed_particles_up=[]
    crossed_particles_down=[]
    rates=open(filename+'.txt','w')
    cumulative_counts=0
    total_residence_time=0
    weighted_residence_time=0
    if use_cutoff==False:
        time_cutoff=total_time
    for ts in universe.trajectory:
        #if ts.frame*timestep<=time_cutoff:    
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
                    #tracking_list[0][water_index]=1
                    tracking_list[water_index]=1
                    #if tracking_list_last[0][water_index]==-2:
                    if tracking_list_last[water_index]==-2:
                        counts_up+=1
                      #  inst_counts_up+=1
                        #pore_time=u.trajectory.frame*timestep - tracking_list[1][water_index]
                     #   total_residence_time+=pore_time
                    #    weighted_residence_time+=permeating_particles*pore_time
                   #     crossed_particles_up.append(water_index)
                elif water_z[water_index]<in_top and water_z[water_index]>in_bottom:
                    #if tracking_list_last[0][water_index]==1 or tracking_list_last[0][water_index]==2:
                    if tracking_list_last[water_index]==1 or tracking_list_last[water_index]==2:
                        #tracking_list[0][water_index]=2
                        tracking_list[water_index]=2
                        #tracking_list[1][water_index]=u.trajectory.frame*timestep
                #    elif tracking_list_last[0][water_index]==-1 or tracking_list_last[0][water_index]==-2:
                    elif tracking_list_last[water_index]==-1 or tracking_list_last[water_index]==-2:
                        tracking_list[water_index]=-2
                        #tracking_list[0][water_index]=-2
                        #tracking_list[1][water_index]=u.trajectory.frame*timestep
                elif water_z[water_index]<in_bottom:
                    #tracking_list[0][water_index]=-1
                    tracking_list[water_index]=-1
                    #if tracking_list_last[0][water_index]==2:
                    if tracking_list_last[water_index]==2:
                        counts_down+=1
                     #   inst_counts_down+=1
                        #pore_time=u.trajectory.frame*timestep - tracking_list[1][water_index]
                        #total_residence_time+=pore_time
                        #weighted_residence_time+=permeating_particles*pore_time
                        #crossed_particles_down.append(water_index)
                #permeating_particles=counts_down+counts_up
                #cumulative_counts+=permeating_particles
                #average_permeation_time=np.divide(total_residence_time,permeating_particles)
                #cumulative_average_permeation_time=np.divide(weighted_residence_time,cumulative_counts)
            if universe.trajectory.frame!=0:
                #inst_rate_up=np.divide(inst_counts_down,u.trajectory.frame*timestep,dtype=np.float64)
                #inst_rate_up
                rate_up=np.divide(counts_up,u.trajectory.frame*timestep,dtype=np.float64)
                #inst_rate_down=np.divide(inst_counts_down,u.trajectory.frame*timestep,dtype=np.float64)
                rate_down=np.divide(counts_down,u.trajectory.frame*timestep,dtype=np.float64)
                rate=np.divide(counts_up-counts_down,u.trajectory.frame*timestep,dtype=np.float64)
                if time_elapsed>=1000:
                    print(rate_down)
                    time_elapsed=0
                inst_rate=np.divide(inst_counts_up-inst_counts_down,u.trajectory.frame*timestep,dtype=np.float64)
                #rates.write(str(rate_up-rate_down)+' '+str(rate_up)+' '+str(rate_down)+' '+str(rate)+ ' '+ str(average_permeation_time) + ' ' + str(cumulative_average_permeation_time) + '\n')
                #rates.write('\n')
                #rates.write(str(inst_rate_up-inst_rate_down)+' '+' '+str(inst_rate_up)+' '+str(inst_rate_down)+' '+str(inst_rate)+'\n')
                #rates.write('\n')
#        else:
    return (rate_up,rate_down,crossed_particles_up,crossed_particles_down)

