import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc
np.set_printoptions(linewidth=5000)    

rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern Roman'],'size':22})
rc('text', usetex=True)

colors = ['#D1655B','#44AA66','#FACD85','#70B9B0','#72B0D7','#E79C5D',
    '#4D5C75','#FFF056','#558C89','#F5CCBA','#A2AB58','#7E8F7C','#005A31']


def read_results(filename):
    lines = []
    with open(filename) as f:
        lines = f.readlines()

    res = []
    for counter, line in enumerate(lines):
        words = line.split()
        if words != []:
            res.append(float(words[0]))

    rows, cols = 6, 4
    res = np.array(res).reshape(cols,rows).T
    # print res
    return res 



def quadrature_benchmark_tabulate(dim=2,compiler="icc"):
    """Tabulate results for latex
    """

    if dim!=2 and dim!=3:
        raise ValueError('Dimension should be either 2 or 3')
    if compiler!="icc" and compiler!="clang" and compiler!="gcc":
        raise ValueError("Compiler should be either 'icc', 'clang' or 'gcc'")
    pwd = os.path.dirname(os.path.realpath(__file__))

    if dim==2:
        res_mn = read_results(pwd+"/2d_benchmark_mooney_rivlin/results_"+compiler)
        res_e1 = read_results(pwd+"/2d_benchmark_electro1/results_"+compiler)
        res_e2 = read_results(pwd+"/2d_benchmark_electro2/results_"+compiler)
    elif dim==3:
        res_mn = read_results(pwd+"/3d_benchmark_mooney_rivlin/results_"+compiler)
        res_e1 = read_results(pwd+"/3d_benchmark_electro1/results_"+compiler)
        res_e2 = read_results(pwd+"/3d_benchmark_electro2/results_"+compiler)

    mm = 6
    for i in range(mm):
        c1 = np.round(res_mn[i,0]/res_mn[i,-1],3)
        c2 = np.round(res_mn[i,1]/res_mn[i,-1],3)
        c3 = np.round(res_mn[i,2]/res_mn[i,-1],3)

        c4 = np.round(res_e1[i,0]/res_e1[i,-1],3)
        c5 = np.round(res_e1[i,1]/res_e1[i,-1],3)
        c6 = np.round(res_e1[i,2]/res_e1[i,-1],3)

        c7 = np.round(res_e2[i,0]/res_e2[i,-1],3)
        c8 = np.round(res_e2[i,1]/res_e2[i,-1],3)
        c9 = np.round(res_e2[i,2]/res_e2[i,-1],3)

        if dim==3:
            if i!=mm-1:
                print "$p="+str(i+1)+"$", "&", c1, "&", c2, "&", c3, "&", c4, "&", c5, "&", c6, "&", c7, "&", c8, "&", c9, r"\\\hline"
            else:
                print "$p="+str(i+1)+"$", "&", c1, "&", c2, "&", c3, "&", c4, "&", c5, "&", c6, "&", c7, "&", c8, "&", c9
        elif dim==2:
            if i!=mm-1:
                print "$p="+str(i+1)+"$", "&", c1, "&", c2, "&", c4, "&", c5, "&", c7, "&", c8, r"\\\hline"
            else:
                print "$p="+str(i+1)+"$", "&", c1, "&", c2, "&", c4, "&", c5, "&", c7, "&", c8
                
    print 


def quadrature_benchmark_plot(dim=2, compiler="icc", model="mn", save=True,figure_name=None):
    """Tabulate results for latex
    """

    if dim!=2 and dim!=3:
        raise ValueError('Dimension should be either 2 or 3')
    if compiler!="icc" and compiler!="clang" and compiler!="gcc":
        raise ValueError("Compiler should be either 'icc', 'clang' or 'gcc'")
    pwd = os.path.dirname(os.path.realpath(__file__))

    if dim==2:
        res_mn = read_results(pwd+"/2d_benchmark_mooney_rivlin/results_"+compiler)
        res_e1 = read_results(pwd+"/2d_benchmark_electro1/results_"+compiler)
        res_e2 = read_results(pwd+"/2d_benchmark_electro2/results_"+compiler)
    elif dim==3:
        res_mn = read_results(pwd+"/3d_benchmark_mooney_rivlin/results_"+compiler)
        res_e1 = read_results(pwd+"/3d_benchmark_electro1/results_"+compiler)
        res_e2 = read_results(pwd+"/3d_benchmark_electro2/results_"+compiler)



    if model=="mn":
        arr1 = res_mn[:,0]/res_mn[:,-1]
        arr2 = res_mn[:,1]/res_mn[:,-1]
        arr3 = res_mn[:,2]/res_mn[:,-1]
    elif model=="e1":
        arr1 = res_e1[:,0]/res_e1[:,-1]
        arr2 = res_e1[:,1]/res_e1[:,-1]
        arr3 = res_e1[:,2]/res_e1[:,-1]
    elif model=="e2":
        arr1 = res_e2[:,0]/res_e2[:,-1]
        arr2 = res_e2[:,1]/res_e2[:,-1]
        arr3 = res_e2[:,2]/res_e2[:,-1]

    if dim==2:
        chaining_percentage     = (arr1 - 1) / (arr1 + arr2 - 2) * 100.
        simd_percentage         = (arr2 - 1) / (arr1 + arr2 - 2) * 100.
    elif dim==3:
        chaining_percentage     = (arr1 - 1) / (arr1 + arr2 + arr3 - 3) * 100
        simd_percentage         = (arr2 - 1) / (arr1 + arr2 + arr3 - 3) * 100
        cross_percentage        = (arr3 - 1) / (arr1 + arr2 + arr3 - 3) * 100

        print chaining_percentage
        print cross_percentage
        print simd_percentage


    mm = 6
    ind = np.arange(mm)
    width = 1.0
    figure, ax = plt.subplots()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])

    if dim==2:
        rects_chain = ax.bar(ind,chaining_percentage+simd_percentage, width, color=colors[8])
        rects_simd = ax.bar(ind,simd_percentage, width, color=colors[5])

        ax.legend((rects_simd[0],rects_chain[0]), 
        ('Vectorisation', 'Operator Chaining'),fontsize=20,loc="upper center", ncol=2, bbox_to_anchor=(0.5, -0.07))


    elif dim==3:
        rects_cross = ax.bar(ind,chaining_percentage+simd_percentage+cross_percentage, width, color=colors[0])
        rects_chain = ax.bar(ind,simd_percentage+chaining_percentage, width, color=colors[8])
        rects_simd  = ax.bar(ind,simd_percentage, width, color=colors[5])

        ax.legend((rects_simd[0],rects_chain[0],rects_cross[0]), 
            ('Vectorisation', 'Operator Chaining','Tensor Cross'), 
            fontsize=20,loc="upper center", ncol=3, bbox_to_anchor=(0.5, -0.07),
            fancybox=True)



    ax.set_xticks((0.45,1.45,2.45,3.45,4.45,5.45))
    ax.set_xticklabels((r"$p=1$",r"$p=2$",r"$p=3$",r"$p=4$",r"$p=5$",r"$p=6$"), fontsize=20)
    plt.ylabel(r"Percentage Contribution",fontsize=24)
    plt.ylim([0,100])

    if save:
        plt.savefig(figure_name,format='eps',dpi=500,bbox_inches='tight',pad_inches=0.1)


    # plt.show()
    plt.close()



if __name__ == "__main__":

    dim=3
    compiler="icc"
    # compiler="clang"
    # compiler="gcc"

    model = "mn"
    # model = "e1"
    # model = "e2"

    # quadrature_benchmark_tabulate(dim=dim,compiler=compiler)


    # quadrature_benchmark_plot(dim=dim,compiler=compiler,model=model)


    for model in ["mn","e1","e2"]:
        for dim in [2,3]:
            folder = "/home/roman/Dropbox/2016_SIMD_Paper/figures/Benchmarks_FEM/"
            figure_name = folder+"/FEM_Bench_"+model+"_"+str(dim)+"D.eps"
            quadrature_benchmark_plot(dim=dim,compiler=compiler,model=model, save=True,figure_name=figure_name)

























# import numpy as np

# def tabulate():

#     filename = "/home/roman/Dropbox/Fastor/benchmark/benchmark_academic/kernel_quadrature/all_results2"
#     with open(filename) as f:
#         lines = f.readlines()

#     each_bench = []
#     for counter, line in enumerate(lines):
#         # print line
#         if "cd" in line:
#             each_bench.append(counter)

#     models = ["2MR", "2MRDS", "2El1", "2El1DS", "2El2", "2El2DS",
#                 "3MR", "3MRDS", "3El1", "3El1DS", "3El2", "3El2DS"]
#     all_speed_ups = {}
#     for i in range(len(each_bench)-0):
#         counter = 0
#         timer = []
#         for j in range(each_bench[i]+2,each_bench[i]+27):
#             if counter % 2 != 0:
#                 lists = lines[j].split(" ")
#                 if lists[7]!='ms.' and lists[7]!='s.':
#                     lists[6] = 10.**(-6)*float(lists[6]) 
#                 elif lists[7]=='ms.':
#                     lists[6] = 10.**(-3)*float(lists[6]) 
#                 if lists[7]=='s.':
#                     lists[6] = float(lists[6])  
#                 timer.append(lists[6])
#             counter+=1
#         timer = np.array(timer)
#         speedup = timer[::2]/timer[1::2]
#         all_speed_ups[models[i]] = speedup

#     for p in range(1,7):
#         print "$p=$"+str(p), "&", np.around(all_speed_ups["2MR"][p-1], decimals=3), "&", \
#                 np.around(all_speed_ups["2El1"][p-1], decimals=3), "&", \
#                 np.around(all_speed_ups["2El2"][p-1], decimals=3), "\\\ \hline"
#     print 
#     for p in range(1,7):
#         print "$p=$"+str(p), "&", np.around(all_speed_ups["3MR"][p-1], decimals=3), "&", \
#                 np.around(all_speed_ups["3El1"][p-1], decimals=3), "&", \
#                 np.around(all_speed_ups["3El2"][p-1], decimals=3), "\\\ \hline"
#     print 
#     for p in range(1,7):
#         print "$p=$"+str(p), "&", np.around(all_speed_ups["3MRDS"][p-1], decimals=3), "&", \
#                 np.around(all_speed_ups["3El1DS"][p-1], decimals=3), "&", \
#                 np.around(all_speed_ups["3El2DS"][p-1], decimals=3), "\\\ \hline"


# if __name__ == "__main__":
#     tabulate()