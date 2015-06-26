import random
import numpy as np
import copy
import MDAnalysis

def rate_sum_update(hopgraph):

    h=hopgraph
    h.filter(exclude={'bulk':True,"outliers":True,"unconnected":True,'oneway':True})
    g=copy.deepcopy(h)
    for node in h.filtered_graph:
        rate_sum=0
        for neighbor in g.filtered_graph[node]:
            rate_sum+=g.filtered_graph[node][neighbor]['k']
        g.filtered_graph[node]['rate_sum']=rate_sum
    rate_sum_bulk=0
    for neighbor in h.graph[1]:
        rate_sum_bulk+=h.graph[1][neighbor]['k']
    g.graph[1]['rate_sum']=rate_sum_bulk
    return g

def flux_calculator(hopgraph,topology,cutoff=1,delta=0.1,steps=1000,particle_num=400,filename='walker_rates',track_entropy_production=False):
    h=hopgraph
    h.filter(exclude={"bulk":True,"outliers":True,"unconnected":True,'oneway':True})
    h=rate_sum_update(h)
    trajectories=np.zeros((4,particle_num))
    u=MDAnalysis.Universe(topology)
    p=u.selectAtoms('protein')
    z_center=p.centerOfGeometry()[2]

    def generate_waiting_time(site):
        rate_sum=h.filtered_graph[site]['rate_sum'] # normalization factor
        r1=random.random()
        waiting_time=-np.log(r1)/rate_sum # assume expontential jump attempt waiting time
        return waiting_time

    def add_particles(trajectories,hopgraph,sites,rate_sum):
        random.seed()
        h=hopgraph
        r1=random.random()
        probability=0
        probability2=0
        added=False
        added_changed=False
        last_site=0
        random.shuffle(sites)
        for site in sites:
            if last_site!=0:
                probability+=h.graph[1][last_site]['k']/rate_sum
            probability2+=h.graph[1][site]['k']/rate_sum
            if probability<=r1<=probability2:
                shuffle_list=range(len(trajectories[0]))
                random.shuffle(shuffle_list)
                shuffle_list_iter=iter(shuffle_list)
                for state in shuffle_list_iter:
                     if trajectories[0][state]==0 and site not in trajectories[0]:
                         trajectories[0][state]=site
                         added=True
                         added_changed=True
                         break
                if added==True:
                    added=False
                    break
            last_site=site
        return trajectories

    def remove_particles(trajectories,h,sites,rate_sum):
        random.seed()
        r1=random.random()
        probability=0
        probability2=0
        removed=False
        last_site=0
        shuffle_list=range(len(trajectories[0]))
        random.shuffle(shuffle_list)
        random_draw_iter=iter(shuffle_list)
        for site in random_draw_iter:
            if trajectories[0][site]!=0 and trajectories[0][site] in sites:
                    if last_site!=0:
                        probability+=h.graph[trajectories[0][last_site]][1]['k']
                    probability2+=h.graph[trajectories[0][site]][1]['k']/rate_sum
                    if probability<=r1<=probability2:
                        removed=True
                        trajectories[0][site]=0
                        trajectories[1][site]=0
                        trajectories[2][site]=0
                        trajectories[3][site]=0
                        break
        last_site=site
        if removed==True:
            removed=False
            return trajectories
        else:
            return False
    def bulk_time_calculator(sites):

        h=hopgraph
        bulk_rate=0
        for site in h.graph[1]:
            if site in sites:
                bulk_rate+=h.graph[1][site]['k']

        bulk_time=1/bulk_rate
        return bulk_time

    def poisson_step(site,site_index,walker_time,waiting_time,delta,entropy_production):
        random.seed()
        w=waiting_time
        r2=random.random()
        probability=0
        for neighbor in h.filtered_graph[site]:
            if type(neighbor)!=str: #filter out rate_sum to prevent bugs
                rate_sum=h.filtered_graph[site]['rate_sum']
                probability+=h.filtered_graph[site][neighbor]['k']
                if probability/rate_sum >= r2:
                    if walker_time>=w:
                            walker_time=0
                            w=generate_waiting_time(neighbor)
                            trajectories[0][site_index]=neighbor
                            trajectories[1][site_index]=0
                            trajectories[2][site_index]=generate_waiting_time(site)
                            if track_entropy_production==True:
                                try:
                                    trajectories[3][site]+=np.log(np.divide(h.filtered_graph[site][neighbor]['k'],h.filtered_graph[neighbor][site]['k']))
                                except:
                                    trajectories[3][site]=np.inf
                            break
                    else:
                        walker_time+=delta
                        trajectories[0][site_index]=site
                        trajectories[1][site_index]=walker_time+delta
                        trajectories[2][site_index]=waiting_time
                        break
        return (trajectories)

    def main(particle_num=particle_num,trajectories=trajectories,up_flux=True,filename=filename):
        random.seed()
        rate_sum_bottom=0
        rate_sum_top=0
        entropy_production=0
        if up_flux==True:
            top_sites_unfiltered = [site for site in h.filtered_graph.nodes() if
            h.site_properties.center[site][2]>=z_center]
            top_sites=[site for site in top_sites_unfiltered if h.graph.has_edge(1,site)]
            bottom_sites_unfiltered = [site for site in h.filtered_graph.nodes() if
            h.site_properties.center[site][2]<=z_center]
            bottom_sites=[site for site in bottom_sites_unfiltered if h.graph.has_edge(site,1)]
        else:
            bottom_sites_unfiltered = [site for site in h.filtered_graph.nodes() if
            h.site_properties.center[site][2]>=z_center]
            bottom_sites=[site for site in top_sites_unfiltered if h.graph.has_edge(1,site)]
            top_sites_unfiltered = [site for site in h.filtered_graph.nodes() if
            h.site_properties.center[site][2]<=z_center]
            top_sites=[site for site in bottom_sites_unfiltered if h.graph.has_edge(site,1)]
        for site in bottom_sites:
                rate_sum_bottom+=h.graph[site][1]['k']
        for site in top_sites:
                rate_sum_top+=h.graph[1][site]['k']

        bulk_time_top=bulk_time_calculator(top_sites)

        bulk_time_bottom=bulk_time_calculator(bottom_sites)
        
        cycle_time_top=bulk_time_top/delta

        cycle_time_bottom=bulk_time_bottom/delta

        steps_elapsed_top=0

        steps_elapsed_bottom=0

        counts=0
        rate=0
        rate_observation=0
        steps_elapsed_total=0
        rates=open(filename+'.txt','w')
        entropy_production=0
        for step in xrange(steps):
            steps_elapsed_total+=1
            random.seed()
            occupied=[site for site in trajectories[0] if site!=0]
            occupied_bottom=[site for site in occupied if site in bottom_sites]
            if steps_elapsed_top>cycle_time_top:
                trajectories=add_particles(trajectories,h,top_sites,rate_sum_top)
                steps_elapsed_top=0
            else:
                steps_elapsed_top+=1
            if steps_elapsed_bottom>cycle_time_bottom:
                out=remove_particles(trajectories,h,bottom_sites,rate_sum_bottom)
                if type(out)!=bool:
                    trajectories=out
                    steps_elapsed_bottom=0
                    counts+=1
            else:
                steps_elapsed_bottom+=1
            shuffle_list=range(len(trajectories[0]))
            random.shuffle(shuffle_list)
            random_draw_iter=iter(shuffle_list)
            for trajectory in random_draw_iter:
                site=trajectories[0][trajectory]
                current_time=trajectories[1][trajectory]
                waiting_time=trajectories[2][trajectory]
                if site!=0.0:
                    trajectories=poisson_step(site,trajectory,current_time,waiting_time,delta,entropy_production)
                    entropy_production+=trajectories[3][site]
            rate=np.divide(counts,delta*step,dtype=np.float64)
            rates.write(str(rate))
            if steps_elapsed_total*delta>=1:
                steps_elapsed_total=0
        rate=np.divide(counts,delta*steps,dtype=np.float64)
        net_entropy_production=np.divide(entropy_production,steps*delta)
        return (trajectories,rate)

    return main(particle_num=particle_num,trajectories=trajectories,filename=filename)

