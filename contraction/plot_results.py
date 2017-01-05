from sys import exit
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc

rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern Roman'],'size':22})
rc('text', usetex=True)

colors = ['#D1655B','#44AA66','#FACD85','#70B9B0','#72B0D7','#E79C5D',
    '#4D5C75','#FFF056','#558C89','#F5CCBA','#A2AB58','#7E8F7C','#005A31']

def read_results(filename):
    lines = []
    with open(filename) as f:
        lines = f.readlines()
    
    arr = []
    for counter, line in enumerate(lines):
        # if counter % 2 is not 0:
        arr.append(float(line))

    return np.array(arr)

def nonisomorphic_tensor_products_plot(no_contracted_indices=1, avx=False, save=False, figure_name=None):

    filename = "products_results"
    scalar_res_gcc = read_results("Scalar_"+filename+"_gcc")
    simd_res_gcc = read_results("SIMD_"+filename+"_gcc")

    scalar_res_clang = read_results("Scalar_"+filename+"_clang")
    simd_res_clang = read_results("SIMD_"+filename+"_clang")

    scalar_res_icc = read_results("Scalar_"+filename+"_icc")
    simd_res_icc = read_results("SIMD_"+filename+"_icc")

    # dim here stands for number of contracting indices
    dim = no_contracted_indices

    if dim==1:
        nn = 1.; stride = np.array([0,1,2,3]) 
    elif dim==2:
        nn = 1.; stride = np.array([4,5,6,7]) 
    elif dim==3:
        nn = 1.; stride = np.array([8,9,10,11]) 
    elif dim==8:
        nn = 1.; stride = np.array([12,13,14,15]) 

    scalar_res_gcc = scalar_res_gcc[stride]*nn
    scalar_res_clang = scalar_res_clang[stride]*nn
    scalar_res_icc = scalar_res_icc[stride]*nn
    simd_res_gcc = simd_res_gcc[stride]
    simd_res_clang = simd_res_clang[stride]
    simd_res_icc = simd_res_icc[stride]

    if not avx:
        scalar_res_gcc = scalar_res_gcc[:2]
        scalar_res_clang = scalar_res_clang[:2]
        scalar_res_icc = scalar_res_icc[:2]
        simd_res_gcc = simd_res_gcc[:2]
        simd_res_clang = simd_res_clang[:2]
        simd_res_icc = simd_res_icc[:2]
    else:
        scalar_res_gcc = scalar_res_gcc[2::]
        scalar_res_clang = scalar_res_clang[2::]
        scalar_res_icc = scalar_res_icc[2::]
        simd_res_gcc = simd_res_gcc[2::]
        simd_res_clang = simd_res_clang[2::]
        simd_res_icc = simd_res_icc[2::]


    N = scalar_res_icc.shape[0]
    ind = np.arange(N)  
    width = 0.3
    gap = 0.08

    fig, ax = plt.subplots()

    rects1 = ax.bar(ind, scalar_res_gcc/simd_res_gcc, width, color=colors[0])
    rects2 = ax.bar(ind+width, scalar_res_clang/simd_res_clang, width, color=colors[8])
    rects3 = ax.bar(ind+2*width, scalar_res_icc/simd_res_icc, width, color=colors[5])


    ax.legend((rects1[0], rects2[0], rects3[0]), 
        ('GCC 6.2.0', 'Clang 3.9.0','ICC 17.0.1'),fontsize=24,loc="best")
    plt.ylabel(r"Speed-up",fontsize=24)

    # draw optimum lines
    # plt.plot([0, 0.9], [4.0, 4.0], color=colors[6], linestyle='-', linewidth=2)
    # plt.plot([1, 1.9], [8.0, 8.0], color=colors[6], linestyle='-', linewidth=2)
    # plt.plot([2, 2.9], [2.0, 2.0], color=colors[6], linestyle='-', linewidth=2)
    # plt.plot([3, 3.9], [4.0, 4.0], color=colors[6], linestyle='-', linewidth=2)


    if not avx:
        plt.text(0.38,0.2,r"SSE",fontsize=18)
        plt.text(1.4,0.2,r"SSE",fontsize=18)
    else:
        plt.text(0.38,0.2,r"AVX",fontsize=18)
        plt.text(1.4,0.2,r"AVX",fontsize=18)

    # plt.ylim([0,20])

    plt.grid('on')
    if dim==1:
        # 1 index
        ax.set_xticks((0.45,1.45))
        if not avx:
            ax.set_xticklabels((r"($2$x$3$x$4$x$5$x$2$) $\times$ ($2$x$3$x$3$x$4$) SP", 
                r"($2$x$3$x$4$x$5$x$2$) $\times$ ($2$x$3$x$3$x$2$) DP"), fontsize=20, rotation=10)
        else:
            ax.set_xticklabels((r"($2$x$3$x$4$x$5$x$2$) $\times$ ($2$x$3$x$3$x$8$) SP", 
                r"($2$x$3$x$4$x$5$x$2$) $\times$ ($2$x$3$x$3$x$4$) DP"), fontsize=20, rotation=10)
    elif dim==2:
        ax.set_xticks((0.45,1.45))
        if not avx:
            ax.set_xticklabels((r"($3$x$4$x$5$x$8$) $\times$ ($3$x$4$x$4$) SP", 
                r"($3$x$4$x$5$x$8$) $\times$ ($3$x$4$x$2$) DP"), fontsize=20, rotation=12)
        else:
            ax.set_xticklabels((r"($3$x$4$x$5$x$8$) $\times$ ($3$x$4$x$8$) SP", 
                r"($3$x$4$x$5$x$8$) $\times$ ($3$x$4$x$8$) DP"), fontsize=20, rotation=12)
    elif dim==3:
        ax.set_xticks((0.45,1.45))
        if not avx:
            ax.set_xticklabels((r"($3$x$4$x$2$x$8$) $\times$ ($3$x$4$x$2$x$4$) SP", 
                r"($3$x$4$x$2$x$8$) $\times$ ($3$x$4$x$2$x$2$) DP"), fontsize=20, rotation=12)
        else:
            ax.set_xticklabels((r"($3$x$4$x$2$x$8$) $\times$ ($3$x$4$x$2$x$8$) SP", 
                r"($3$x$4$x$2$x$8$) $\times$ ($3$x$4$x$2$x$8$) DP"), fontsize=20, rotation=12)
    elif dim==8:
        ax.set_xticks((0.45,1.55))
        if not avx:
            ax.set_xticklabels((r"($2$x$3$x$2$x$3$x$2$x$3$x$2$x$3$x$2$) $\times$ ($2$x$3$x$2$x$3$x$2$x$3$x$2$x$3$x$4$) SP",
                r"($2$x$3$x$2$x$3$x$2$x$3$x$2$x$3$x$2$) $\times$ ($2$x$3$x$2$x$3$x$2$x$3$x$2$x$3$x$2$) DP"), 
                fontsize=18, rotation=5)
        else:
            ax.set_xticklabels((r"($2$x$3$x$2$x$3$x$2$x$3$x$2$x$3$x$2$) $\times$ ($2$x$3$x$2$x$3$x$2$x$3$x$2$x$3$x$4$) SP",
                r"($2$x$3$x$2$x$3$x$2$x$3$x$2$x$3$x$2$) $\times$ ($2$x$3$x$2$x$3$x$2$x$3$x$2$x$3$x$2$) DP"), 
                fontsize=18, rotation=5)

    if save:
        print(figure_name)
        plt.savefig(figure_name,format='eps',dpi=500,bbox_inches='tight',pad_inches=0.1)

    plt.show()


if __name__ == "__main__":
    
    no_contracted_indices = 1
    avx = True

    folder = "/home/roman/Dropbox/2016_SIMD_Paper/figures/Benchmarks_SIMD/"
    if not avx:
        figure_name = folder+"Benchmark_SSE_Contraction_Index_{}.eps".format(no_contracted_indices)
    else:
        figure_name = folder+"Benchmark_AVX_Contraction_Index_{}.eps".format(no_contracted_indices)

    nonisomorphic_tensor_products_plot(no_contracted_indices=no_contracted_indices, avx=False)
    # nonisomorphic_tensor_products_plot(no_contracted_indices=no_contracted_indices, avx=avx, save=True, figure_name=figure_name)