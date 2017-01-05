import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator
import matplotlib.cm as cm
from matplotlib import rc
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern Roman'],'size':20})
rc('text', usetex=True)

colors = ['#D1655B','#44AA66','#FACD85','#70B9B0','#72B0D7','#E79C5D',
    '#4D5C75','#FFF056','#558C89','#F5CCBA','#A2AB58','#7E8F7C','#005A31']



def to_eng_str(x):
    y = np.abs(x)
    exponent = int(np.log10(y))
    return '$'+str(np.round(y/10.**exponent,3))+r'$$\times 10^{'+str(exponent)+r'}$'


def read_results(filename):
    lines = []
    with open(filename) as f:
        lines = f.readlines()

    flop_nodp = []; flop_dp =[]; temp_size = []
    cost = []
    for counter, line in enumerate(lines):
        words = line.split()
        if "FLOP" in line:
            flop_nodp.append(int(words[5]))
            flop_dp.append(int(words[10]))
            temp_size.append(int(words[-1]))
        else:
            if words:
                cost.append(float(words[0]))

    return np.array(flop_nodp), np.array(flop_dp), np.array(temp_size), np.array(cost) 



def operation_minimisation_plots(cache_idx=0, ntemp=1, save=False, figure_name=None):

    from scipy.stats import itemfreq

    if ntemp==1:
        flop_nodp, flop_dp, temp_size, cost_dp = read_results("gcc_res_dp")
        cost_nodp = read_results("gcc_res_nodp")[-1]
    else:
        flop_nodp, flop_dp, temp_size, cost_dp = read_results("gcc_res_2_dp")
        cost_nodp = read_results("gcc_res_2_nodp")[-1]

    flops_saved = flop_nodp - flop_dp

    if ntemp==1:
        temp_size[temp_size==10] = 10*1024
        temp_size[temp_size==20] = 20*1024
        temp_size[temp_size==80] = 80*1024

    tmp_idx = []
    counter = 0
    freqs = itemfreq(temp_size)
    for i in range(freqs.shape[0]):
        tmp_idx.append(np.arange(counter,counter+freqs[i,1]))
        counter = np.sum(freqs[:i+1,1])
    tt = cache_idx



    caches = ["0_5xL1","L1","0_5xL2","L2","0_5xL3","L3","4xL3"]
    cache_labels = ["L1","L1","L2","L2","L3","L3","Memory"] 

    fig, ax = plt.subplots()

    N = tmp_idx[tt].shape[0]
    ind = np.arange(N)  
    width = 0.7
    color = colors[10]

    speed_ups = cost_nodp[tmp_idx[tt]]/cost_dp[tmp_idx[tt]]
    rects = ax.bar(ind+0.1, speed_ups, width, color=color, log=True)   

    reduced_flops = [to_eng_str(number) for number in flops_saved[tmp_idx[tt]]]
    ax.set_xticks(ind+0.3)
    ax.set_xticklabels(reduced_flops, rotation=40, fontsize=18, ha='center')

    plt.ylabel(r"Speed-up (for data in "+cache_labels[tt]+")",fontsize=22)

    plt.grid('on')
    xlim = ax.get_xlim()
    ax.set_xlim([xlim[0]+0.1,xlim[1]])
    ylim = ax.get_ylim()
    plt.ylim([ylim[0],1000])

    h_label = mpatches.Patch(color=color, label=r'Memory/FLOP Optimal Contraction')
    plt.legend(handles=[h_label], fontsize=18)

    for rect, label in zip(rects, speed_ups):
        height = rect.get_height()
        plt.text(rect.get_x() + rect.get_width()/2, height+.0, np.round(label,3), ha='center', va='bottom', fontsize=18)

    if save:
        if ntemp==1:
            plt.savefig(figure_name+caches[tt]+".eps",format='eps',dpi=500,bbox_inches='tight',pad_inches=0.1)
        elif ntemp==2:
            plt.savefig(figure_name+caches[tt]+".eps",format='eps',dpi=500,bbox_inches='tight',pad_inches=0.1)

    # plt.show()
    plt.close()



if __name__ == "__main__":

    folder = "/home/roman/Dropbox/2016_SIMD_Paper/figures/Benchmarks_DepthFirst/"
    ntemp=1
    if ntemp==1:
        figure_name = folder+"Benchmark_DepthFirst_OneTemp_"
    elif ntemp==2:
        figure_name = folder+"Benchmark_DepthFirst_TwoTemps_"

    for i in range(7):     
        # operation_minimisation_plots(cache_idx=i,ntemp=ntemp)
        operation_minimisation_plots(cache_idx=i, ntemp=ntemp, save=True, figure_name=figure_name)
















































# import numpy as np
# import matplotlib as mpl
# import matplotlib.pyplot as plt
# import matplotlib.cm as cm
# from matplotlib import rc

# rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern Roman'],'size':20})
# rc('text', usetex=True)
# # rc('axes',color_cycle=['#D1655B','#FACD85','#72B0D7','#E79C5D','#4D5C75','#E79C5D'])

# colors = ['#D1655B','#44AA66','#FACD85','#70B9B0','#72B0D7','#E79C5D',
#     '#4D5C75','#FFF056','#558C89','#F5CCBA','#A2AB58','#7E8F7C','#005A31']

# def read_results(filename):
#     lines = []
#     with open(filename) as f:
#         lines = f.readlines()

#     flop_nodp = []; flop_dp =[]; temp_size = []
#     cost_nodp = []; cost_dp = [];
#     for counter, line in enumerate(lines):
#         words = line.split()
#         if "FLOP" in line:
#             flop_nodp.append(int(words[5]))
#             flop_dp.append(int(words[10]))
#             temp_size.append(int(words[-1]))
#         else:
#             if words:
#                 cost_nodp.append(float(words[0]))
#                 cost_dp.append(float(words[1])) 

#     return np.array(flop_nodp), np.array(flop_dp), np.array(temp_size), np.array(cost_nodp), np.array(cost_dp) 



# def operation_minimisation_plots(cache_idx=0, ntemp=1, save=False, figure_name=None):

#     from scipy.stats import itemfreq

#     if ntemp==1:
#         flop_nodp, flop_dp, temp_size, cost_nodp, cost_dp = read_results("gcc_res_dp")
#     else:
#         flop_nodp, flop_dp, temp_size, cost_nodp, cost_dp = read_results("gcc_res_2")

#     flops_saved = flop_nodp - flop_dp

#     if ntemp==1:
#         temp_size[temp_size==10] = 10*1024
#         temp_size[temp_size==20] = 20*1024
#         temp_size[temp_size==80] = 80*1024

#     tmp_idx = []
#     counter = 0
#     freqs = itemfreq(temp_size)
#     for i in range(freqs.shape[0]):
#         tmp_idx.append(np.arange(counter,counter+freqs[i,1]))
#         counter = np.sum(freqs[:i+1,1])
#     # fix
#     # if ntemp==2:
#     #     tmp_idx[0] = tmp_idx[0][1:]
#     #     tmp_idx[1] = tmp_idx[1][1:]
#     #     # print tmp_idx[0]

#     tt = cache_idx
#     # tt corresponds to 

#     # if ntemp==2:
#     #     if tt==0:
#     #         cost_dp[1] /=2.5 


#     caches = ["0_5xL1","L1","0_5xL2","L2","0_5xL3","L3","4xL3"] 

#     fig, ax = plt.subplots()

#     N = tmp_idx[tt].shape[0]
#     ind = np.arange(N)  
#     width = 0.3
#     gap = 0.1
#     rects1 = ax.bar(ind, cost_dp[tmp_idx[tt]], width, color=colors[8], log=True)
#     rects2 = ax.bar(ind+width, cost_nodp[tmp_idx[tt]], width, color=colors[5], log=True)
        
#     # exit()

#     ax.set_xticks(ind+0.1)
#     ax.set_xticklabels((flops_saved[tmp_idx[tt]]), rotation=40, fontsize=18, ha='center')
#     plt.ylim([10**(-6),100])

#     plt.ylabel(r"Run-Time Per Call (sec)",fontsize=22)

#     plt.grid('on')

#     ax.legend((rects1[0], rects2[0]), 
#         (r'FLOP Optimal Contraction',
#             r'Memory Optimal Contraction'),fontsize=20,loc="best")

#     if save:
#         if ntemp==1:
#             plt.savefig(folder+"Benchmark_DepthFirst_OneTemp_"+caches[tt]+".eps",format='eps',dpi=500,bbox_inches='tight',pad_inches=0.1)
#         elif ntemp==2:
#             plt.savefig(folder+"Benchmark_DepthFirst_TwoTemps_"+caches[tt]+".eps",format='eps',dpi=500,bbox_inches='tight',pad_inches=0.1)

#     plt.show()



# if __name__ == "__main__":
#     # operation_minimisation_plots()

#     folder = "/home/roman/Dropbox/2016_SIMD_Paper/figures/Benchmarks_DepthFirst/"

#     for i in range(7):
                

#         # do_plot(cache_idx=i)
#         # do_plot(cache_idx=i,ntemp=2,save=True)

#         operation_minimisation_plots(cache_idx=i,ntemp=1,save=False)