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
        if counter % 2 is not 0:
            arr.append(float(line))

    return np.array(arr)


def isomorphic_tensor_products_plot(tensor_dimension=2, save=False, figure_name=None):
    """Plots results of isomorphic tensor products 
    """

    filename = "products_results"
    scalar_res_gcc = read_results("Scalar_"+filename+"_gcc")
    simd_res_gcc = read_results("SIMD_"+filename+"_gcc")

    scalar_res_clang = read_results("Scalar_"+filename+"_clang")
    simd_res_clang = read_results("SIMD_"+filename+"_clang")

    scalar_res_icc = read_results("Scalar_"+filename+"_icc")
    simd_res_icc = read_results("SIMD_"+filename+"_icc")

    dim = tensor_dimension

    if dim==2:
        stride = np.array([0,1,2,3]) # 2D
    elif dim==3:
        stride = np.array([4,5,6,7]) # 3D
    elif dim==4:
        stride = np.array([8,9,10,11]) # 4D
    elif dim==5:
        stride = np.array([12,13,14,15]) # 5D
    elif dim==6:
        stride = np.array([16,17,18,19]) # 6D
    elif dim==8:
        stride = np.array([20,21,22,23]) # 8D

    scalar_res_gcc = scalar_res_gcc[stride]
    scalar_res_clang = scalar_res_clang[stride]
    scalar_res_icc = scalar_res_icc[stride]
    simd_res_gcc = simd_res_gcc[stride]
    simd_res_clang = simd_res_clang[stride]
    simd_res_icc = simd_res_icc[stride]

    N = scalar_res_icc.shape[0]
    ind = np.arange(N)  
    width = 0.3
    gap = 0.08

    fig, ax = plt.subplots()

    to_plot_gcc = scalar_res_gcc/simd_res_gcc
    to_plot_clang = scalar_res_clang/simd_res_clang
    to_plot_icc = scalar_res_icc/simd_res_icc

    # rects1 = ax.bar(ind, to_plot_gcc, width, color=colors[2])
    # rects2 = ax.bar(ind+width, to_plot_clang, width, color=colors[0])
    # rects3 = ax.bar(ind+2*width, to_plot_icc, width, color=colors[4])

    rects1 = ax.bar(ind, to_plot_gcc, width, color=colors[0])
    rects2 = ax.bar(ind+width, to_plot_clang, width, color=colors[8])
    rects3 = ax.bar(ind+2*width, to_plot_icc, width, color=colors[5])


    ax.legend((rects1[0], rects2[0], rects3[0]), 
        ('GCC 6.2.0', 'Clang 3.9.0','ICC 17.0.1'),fontsize=24)
    # plt.ylabel(r"Speed-up $f(FLOP/cycle,\;byte/cycle)$",fontsize=18)
    plt.ylabel(r"Speed-up",fontsize=24)

    # draw optimum lines
    # plt.plot([0, 0.9], [4.0, 4.0], color=colors[6], linestyle='-', linewidth=2)
    # plt.plot([1, 1.9], [8.0, 8.0], color=colors[6], linestyle='-', linewidth=2)
    # plt.plot([2, 2.9], [2.0, 2.0], color=colors[6], linestyle='-', linewidth=2)
    # plt.plot([3, 3.9], [4.0, 4.0], color=colors[6], linestyle='-', linewidth=2)


    plt.text(0.315,0.2,r"SSE",fontsize=20)
    plt.text(1.3,0.2,r"AVX",fontsize=20)
    plt.text(2.314,0.2,r"SSE",fontsize=20)
    plt.text(3.3,0.2,r"AVX",fontsize=20)

    # plt.ylim([0,20])

    plt.grid('on')
    if dim==2:
        # 2D
        ax.set_xticks((0.48,1.45,2.5,3.5))
        ax.set_xticklabels((r"$4$x$4$ SP", r"$2$x$16$ SP", r"$3$x$2$ DP", r"$4$x$4$ DP"), 
            fontsize=22, rotation=30) # 2D
    elif dim==3:
        # 3D
        ax.set_xticks((0.45,1.45,2.5,3.5))
        ax.set_xticklabels((r"$4$x$4$x$4$ SP", r"$2$x$3$x$16$ SP", r"$4$x$3$x$2$ DP", 
            r"$2$x$3$x$4$ DP"), 
            fontsize=22, rotation=28) # 3D
    elif dim==4:
        # 4D
        ax.set_xticks((0.45,1.45,2.5,3.5))
        ax.set_xticklabels((r"$2$x$3$x$4$x$4$ SP", r"$2$x$3$x$4$x$8$ SP", 
            r"$5$x$4$x$3$x$2$ DP", r"$2$x$3$x$5$x$4$ DP"), 
            fontsize=20, rotation=26) 
    elif dim==5:
        # 5D
        ax.set_xticks((0.45,1.45,2.5,3.5))
        ax.set_xticklabels((r"$2$x$2$x$2$x$2$x$4$ SP", r"$2$x$2$x$2$x$3$x$16$ SP", 
            r"$2$x$2$x$2$x$2$x$2$ DP", r"$2$x$2$x$2$x$3$x$16$ DP"), 
            fontsize=20, rotation=24) 
    elif dim==6:
        # 6D
        ax.set_xticks((0.45,1.45,2.5,3.5))
        ax.set_xticklabels((r"$2$x$2$x$2$x$2$x$2$x$4$ SP", r"$2$x$2$x$2$x$2$x$2$x$8$ SP",
            r"$2$x$2$x$2$x$2$x$2$x$2$ DP", r"$2$x$2$x$2$x$2$x$2$x$4$ DP"), 
            fontsize=20, rotation=20) 
    elif dim==8:
        # 8D
        ax.set_xticks((0.45,1.45,2.5,3.5))
        ax.set_xticklabels((r"$2$x$2$x$2$x$2$x$2$x$2$x$2$x$4$ SP",
            r"$2$x$2$x$2$x$2$x$2$x$2$x$8$ SP", r"$2$x$2$x$2$x$2$x$2$x$2$x$2$ DP", 
            r"$2$x$2$x$2$x$2$x$2$x$2$x$4$ DP"), 
            fontsize=20, rotation=16) 

    if save:
        plt.savefig(figure_name,format='eps',dpi=500,bbox_inches='tight',pad_inches=0.1)

    plt.show()


if __name__ == "__main__":
    tensor_dimension = 8
    folder = "/home/roman/Dropbox/2016_SIMD_Paper/figures/Benchmarks_SIMD/"
    figure_name = folder+"Benchmark_SIMD_{}D.eps".format(tensor_dimension)
    isomorphic_tensor_products_plot(tensor_dimension=tensor_dimension, save=False)
    # isomorphic_tensor_products_plot(tensor_dimension=tensor_dimension, save=True, figure_name=figure_name)